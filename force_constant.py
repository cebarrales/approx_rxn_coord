# -*- coding: utf-8 -*-
# PyGopt: Python Geometry Optimization.
# Copyright (C) 2011-2018 The HORTON/PyGopt Development Team
#
# This file is part of PyGopt.
#
# PyGopt is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# PyGopt is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
""""Class for parsing data from FCHK file"""

import numpy as np
from numpy import matrix
import sys


class FCHKFile(object):
    """Reader for Formatted checkpoint files

       After initialization, the data from the file is available in the fields
       dictionary. Also the following attributes are read from the file: title,
       command, lot (level of theory) and basis.
    """

    def __init__(self, filename: str, ignore_errors: bool = False, field_labels=None):
        """
           Arguments:
            | ``filename``  --  The formatted checkpoint file

           Optional arguments:
            | ``ignore_errors``  --  Try to read incorrectly formatted files
                                     without raising exceptions [default=False]
            | ``field_labels``  --  When given, only these fields are read from
                                    the formatted checkpoint file. (This can
                                    save a lot of time.)
        """
        self.filename = filename
        self.ignore_errors = ignore_errors
        try:
            if field_labels is not None:
                field_labels = set(field_labels)
                field_labels.add("Atomic numbers")
                field_labels.add("Current cartesian coordinates")
            self._read(filename, field_labels)
        except FileFormatError:
            if ignore_errors:
                pass
            else:
                raise
        self._analyze()

    def _read(self, filename, field_labels=None):
        """Read all the requested fields

           Arguments:
            | ``filename``  --  the filename of the FCHK file
            | ``field_labels``  --  when given, only these fields are read
        """

        # if fields is None, all fields are read
        def read_field(file):
            """Read a single field"""
            datatype = None
            while datatype is None:
                # find a sane header line
                line = file.readline()
                if line == "":
                    return False

                label = line[:43].strip()
                if field_labels is not None:
                    if len(field_labels) == 0:
                        return False
                    elif label not in field_labels:
                        return True
                    else:
                        field_labels.discard(label)
                line = line[43:]
                words = line.split()
                if len(words) == 0:
                    return True

                if words[0] == 'I':
                    datatype = int
                    unreadable = 0
                elif words[0] == 'R':
                    datatype = float
                    unreadable = np.nan

            if len(words) == 2:
                try:
                    value = datatype(words[1])
                except ValueError:
                    return True
            elif len(words) == 3:
                if words[1] != "N=":
                    raise FileFormatError(
                        "Unexpected line in formatted checkpoint file %s\n%s" %
                        (filename, line[:-1]))
                length = int(words[2])
                value = np.zeros(length, datatype)
                counter = 0
                try:
                    while counter < length:
                        line = file.readline()
                        if line == "":
                            raise FileFormatError(
                                "Unexpected end of formatted checkpoint file %s"
                                % filename)
                        for word in line.split():
                            try:
                                value[counter] = datatype(word)
                            except (ValueError, OverflowError):
                                print('WARNING: could not interpret word while reading %s: %s' % (
                                    word, self.filename))
                                if self.ignore_errors:
                                    value[counter] = unreadable
                                else:
                                    raise
                            counter += 1
                except ValueError:
                    return True
            else:
                raise FileFormatError(
                    "Unexpected line in formatted checkpoint file %s\n%s" %
                    (filename, line[:-1]))

            self.fields[label] = value
            return True

        self.fields = {}
        with open(filename, 'r') as f:
            self.title = f.readline()[:-1].strip()
            words = f.readline().split()
            if len(words) == 3:
                self.command, self.lot, self.basis = words
            elif len(words) == 2:
                self.command, self.lot = words
            else:
                raise FileFormatError(
                    'The second line of the FCHK file should contain two or three words.'
                )

            while read_field(f):
                pass

            f.close()

    def _analyze(self):
        """Convert a few elementary fields into a molecule object"""
        if ("Atomic numbers" in self.fields) and (
                "Current cartesian coordinates" in self.fields):
            self.molecule = Molecule(
                self.fields["Atomic numbers"],
                np.reshape(self.fields["Current cartesian coordinates"],
                           (-1, 3)),
                self.title, )

    def get_coordinates(self):
        """Return cartesian coordinates in format check file"""
        return self.fields.get("Current cartesian coordinates")

    def get_hessian(self):
        """Return the hessian"""
        force_const = self.fields.get("Cartesian Force Constants")
        if force_const is None:
            return None
        N = len(self.molecule.numbers)
        result = np.zeros((3 * N, 3 * N), float)
        counter = 0
        for row in range(3 * N):
            result[row, :row + 1] = force_const[counter:counter + row + 1]
            result[:row + 1, row] = force_const[counter:counter + row + 1]
            counter += row + 1
        return result

    def get_gradient(self):
        return self.fields.get("Cartesian Gradient")

    def get_energy(self):
        return self.fields.get("Total Energy")


class Molecule(object):
    def __init__(self, numbers, coordinates, title):
        self.numbers = numbers
        self.coordinates = coordinates
        self.title = title


class FileFormatError(Exception):
    pass


def force_constant(r, r_1, p_1, p):
    """ Obtain the force constant of R and P from the Hessian of both and
        the next and previous point of an IRC respectively """
    x_0 = FCHKFile(r)
    x_1 = FCHKFile(r_1)
    x_p_1 = FCHKFile(p_1)
    x_p = FCHKFile(p)
    hessian_r = matrix(x_0.get_hessian())
    hessian_p = matrix(x_p.get_hessian())
    coord_r = matrix(x_0.get_coordinates())
    coord_r_1 = matrix(x_1.get_coordinates())
    coord_p_1 = matrix(x_p_1.get_coordinates())
    coord_p = matrix(x_p.get_coordinates())
    dcoord_r = coord_r_1 - coord_r
    dcoord_p = coord_p - coord_p_1
    k_r = (dcoord_r * hessian_r * dcoord_r.T) / (dcoord_r * dcoord_r.T)
    k_p = (dcoord_p * hessian_p * dcoord_p.T) / (dcoord_p * dcoord_p.T)
    return float(k_r), float(k_p)


def ts_force_constant(ts_file):
    """ Read TS force constant from .log file of a TS calculation """
    f = open(ts_file, "r")
    fc = []
    for line in f:
        if "Frc consts  --" in line:
            fc.append(line.split()[3])
    f_const_ts = float(fc[0]) * -0.22937
    return f_const_ts

