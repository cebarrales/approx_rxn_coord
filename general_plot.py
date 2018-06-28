#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
This module allows to plot the results from the coordinate transformation from IRC
to QRC and SRC. It reads the quartic.dat and sextic.dat of the coordinate transformation.
The final plot will have the QRC and SCR profiles in order to compare them.
The name of new file with the plot must be specified in the first parameter of the function,
if you want to save the picture you must give "yes" as the second parameter and if you want
to include the IRC in order to compare it with the QRC and SRC, you must write "yes" as the
third parameter of the function compare_plot().

Example:

compare_plot('test', 'yes', 'yes')

It will save a picture called test.png in the folder where you execute the code and it will
contain a QRC, SCR and IRC plot.

"""

from sympy import *
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import font_manager
import sys
import re


DEFAULT_NAME = sys.argv[3]
WANT_PLT = 'no_save'
IRC = 'no'
quartic_file = sys.argv[1]
sextic_file = sys.argv[2]


def max_points(quartic_data):
    """ Obtains the number of points of IRC calculation """
    file = open(quartic_data, "r")
    for num, line in enumerate(file, 1):
        if "1.0\t" in line:
            max_point_number = num - 53
            return max_point_number


def energy(quartic_data):
    """ Save the energy on each point of IRC calculation in a list"""
    file = open(quartic_data, "r")
    energies = []
    for line in file:
        if "IRC\t" in line:
            for _ in list(range(2)):
                next(file)
            for _ in list(range(max_points(quartic_data))):
                energy_k = ((next(file)).split()[2])
                energies.append(float(energy_k))
    return energies


def qrc(quartic_data):
    """ Save the energy on each point of IRC calculation in a list"""
    file = open(quartic_data, "r")
    qrc_points = []
    for line in file:
        if "IRC\t" in line:
            for _ in list(range(2)):
                next(file)
            for _ in list(range(max_points(quartic_data))):
                qrc_k = ((next(file)).split()[1])
                qrc_points.append(float(qrc_k))
    return qrc_points


def src(sextic_data, quartic_data):
    """ Save the energy on each point of IRC calculation in a list"""
    file = open(sextic_data, "r")
    src_points = []
    for line in file:
        if "IRC\t" in line:
            for _ in list(range(2)):
                next(file)
            for _ in list(range(max_points(quartic_data))):
                src_k = ((next(file)).split()[1])
                src_points.append(float(src_k))
    return src_points


def irc(quartic_data):
    """ Save the energy on each point of IRC calculation in a list"""
    file = open(quartic_data, "r")
    irc_points = []
    for line in file:
        if "IRC\t" in line:
            for _ in list(range(2)):
                next(file)
            for _ in list(range(max_points(quartic_data))):
                irc_k = ((next(file)).split()[0])
                irc_points.append(float(irc_k))
    return irc_points


def compare_plot(save_plot=WANT_PLT, include_irc=IRC, name=DEFAULT_NAME):
    """ This function will plot the QRC, SRC and IRC (if it is required) """
    e = energy(quartic_file)
    q = qrc(quartic_file)
    s = src(sextic_file, quartic_file)
    i = irc(quartic_file)
    plt.rcParams['font.family'] = 'Times New Roman'
    if include_irc == 'yes':
        plt.plot(i, e, label='IRC', marker='o', fillstyle='none',
                 markersize='2.5', linewidth=0.5, c='g')
        plt.plot(q, e, label='QRC', marker='o', fillstyle='none',
                 markersize='3', linewidth=0.5)
        plt.plot(s, e, label='SRC', marker='o', fillstyle='none',
                 markersize='3', linewidth=0.5)
        plt.title('QRC vs. SRC vs. IRC', fontsize=13)
        plt.ylabel('Potential Energy (kcal/mol)', fontsize=13)
        plt.xlabel('Reaction coordinate', fontsize=13)
        plt.legend(fontsize=12)
    else:
        plt.plot(q, e, label='QRC', marker='o', fillstyle='none',
                 markersize='3', linewidth=0.5)
        plt.plot(s, e, label='SRC', marker='o', fillstyle='none',
                 markersize='3', linewidth=0.5)
        plt.title(name, fontsize=13, fontweight='bold')
        plt.ylabel('Potential Energy (kcal/mol)', fontsize=13)
        plt.xlabel('Reaction coordinate', fontsize=13)
        plt.legend(fontsize=12)
    if save_plot == 'yes':
        return plt.savefig(name + '.png', dpi=1080)
    else:
        plt.show()
    return


compare_plot('yes', 'no')
