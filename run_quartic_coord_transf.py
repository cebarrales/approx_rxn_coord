#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Quartic reaction coordinate
    
    This module contains a set of functions, which allow to perform an optimization of the coefficients
    in the polinomial expansion of quartic reaction coordinate (QRC) model (a*x**4 + b*x**3 + c*x**2),
    taking data from a file containing the reaction coordinate and energy in each point of an IRC calculation
    and performs a coordinate transformation from IRC to QRC.
    
    Example
    -------
    You must have a file with the IRC and energy values of each point. From terminal execute:
    
    run_quartic_coord_transf energy.dat

    The result will be a file called quartic.dat with two sections: The first one is a summary of the main r
    esults from the new QRC and the second one is the new coordinate in each point of IRC.
    
    """

from quartic_reaction_coordiante import *
from sympy import *
import sys

irc_file = sys.argv[1]


def max_points(irc_data):
    """Save the total number of points in IRC """
    file = open(irc_data, "r")
    max_points_irc = (int(file.readlines()[-1].split()[0]))
    return max_points_irc


def energy(irc_data):
    """Save the energy in each points of IRC """
    file = open(irc_data, "r")
    energies = []
    for line in file:
        if "Reaction coord              Energy" in line:
            for _ in list(range(max_points(irc_data))):
                energy_k = (float(next(file).split()[2]))
                energies.append(energy_k)
    return energies


def coord_irc(irc_data):
    print('Analyzing each point of IRC ...')
    """Save the coordinate value of each points of IRC """
    file = open(irc_data, "r")
    coordinate = []
    coor_zero_to_one = []
    for line in file:
        if "Reaction coord              Energy" in line:
            for _ in list(range(max_points(irc_data))):
                coord_k = (float(next(file).split()[1]))
                coordinate.append(coord_k)
        for coord in coordinate:
            a_coord = round((coord - coordinate[0]) /
                            (coordinate[max_points(irc_data) - 1] - coordinate[0]), 4)
            coor_zero_to_one.append(a_coord)
    file.close()
    print('Done!')
    return coor_zero_to_one


def coord_transformation(irc_data):
    """Perform a coordinate transformation from IRC to QRC """
    file = open(irc_data, "r")
    output_file = (open("quartic.dat", "w"))
    eact, erxn = [float(max(energy(irc_data))),
                  float(file.readlines()[-1].split()[2])]
    a, b, c = quarticrxn(eact, erxn, "yes", output_file)
    print('Transforming coordinate from IRC to QRC ...')
    x = Symbol('x')
    sb = Matrix([0.1])
    x_transformed = []
    iteration = 0
    for i in energy(irc_data):
        f = a * x ** 4 + b * x ** 3 + c * x ** 2 - i
        m = Matrix([f])
        j_i = (m.jacobian([x])) ** (-1)
        if i == 0:
            s0 = sb
            while m.subs([(x, s0[0])]).norm() > 5e-8:
                s0 = s0 - (m.subs([(x, s0[0])]) * j_i.subs([(x, s0[0])]))
                iteration = iteration + 1
                # print(s0[0])
            x_transformed.append(round(s0[0], 4))
            sb = Matrix([s0[0] + 1 / (max_points(irc_data) + 1)])
        elif i == erxn:
            while m.subs([(x, sb[0])]).norm() > 5e-10:
                sb = sb - (m.subs([(x, sb[0])]) * j_i.subs([(x, sb[0])]))
                iteration = iteration + 1
                # print(sb[0])
            x_transformed.append(round(sb[0], 4))
            sb = Matrix([0.9999])
        else:
            while m.subs([(x, sb[0])]).norm() > 5e-10:
                sb = sb - (m.subs([(x, sb[0])]) * j_i.subs([(x, sb[0])]))
                iteration = iteration + 1
                # print(sb[0])
            x_transformed.append(round(sb[0], 4))
            if sb[0] > 0.99:
                sb = Matrix([0.999])
            else:
                sb = Matrix([sb[0] + 1 / (max_points(irc_data + 1))])
    print('Done!')
    output_file.close()
    return x_transformed


def save_data(irc_data):
    """Save data of QRC analysis and the coordinate transformation """
    output_file = (open("quartic.dat", "a"))
    irc = coord_irc(irc_data)
    qrc = coord_transformation(irc_data)
    energy_of_each_point = energy(irc_data)
    print('Saving data ...')
    quartic_plot(energy_of_each_point, irc, qrc, "plot", "save")
    output_file.write('\n=====================================================\n\n')
    output_file.write('\n              Coordinate transformation              \n')
    output_file.write('                   from IRC to QRC                   \n\n')
    output_file.write('\n=====================================================\n\n')
    output_file.write("\t IRC\t\t  QRC\t\t  Energy\n")
    output_file.write("\t\t\t  \t\t(kcal/mol)\n")
    output_file.write('=====================================================\n')
    for point in list(range(max_points(irc_file))):
        output_file.write("\t" + str(irc[point]) + "\t\t" + str(qrc[point]) + "\t\t" +
                          str(energy_of_each_point[point]) + "\n")
    print('Done!')
    output_file.write('\n==========================END========================\n')
    output_file.close()


save_data(irc_file)
