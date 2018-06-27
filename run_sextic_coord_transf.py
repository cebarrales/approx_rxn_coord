#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Sextic reaction coordinate
    
    This module contains a set of functions, which allow to perform an optimization of the coefficients
    in the polinomial expansion of quartic reaction coordinate (QRC) model (a*x**6 + b*x**5 + c*x**4 + d*x**3 + e*x**2),
    taking data from a file containing the reaction coordinate and energy in each point of an IRC calculation
    and performs a coordinate transformation from IRC to SRC. It requires the force constant of R and P, therefore
    it is necessary additional files that you must have:
    1) A .log file with the optimization of the TS (with freq calculation).
    2) A .fchk file from a freq calculation of the R and P.
    3) A .fchk file from a single point calculation of the R+1 and P-1 points of IRC.
    4) A file containing the coordinate and energy of each point of IRC.
    
    Example
    -------
    Having all the files required, from terminal execute:
    
    run_sextic_coord_transf energy.dat ts_file.log sp_0.fchk sp_1.fchk sp_p_1.fchk sp_p.fchk

    Results will be saved in the sextic.dat file and the plot will be saved in the plot.png file.
    
    """

from sextic_reaction_coordinate import *
from force_constant import force_constant, ts_force_constant
from sympy import *
import sys


irc_file = sys.argv[1]
ts_file = sys.argv[2]
irc_0 = sys.argv[3]
irc_1 = sys.argv[4]
irc_p_1 = sys.argv[5]
irc_p = sys.argv[6]

def max_points(irc_file):
    """Save the total number of points in IRC """
    file = open(irc_file,"r")
    max_points_irc = (int(file.readlines()[-1].split()[0])) # del archivo, lea la ultima linea (-1) y entregue la primera columna (1) como entero.
    return max_points_irc



def energy(irc_file):
    """Save the energy in each points of IRC """
    file = open(irc_file,"r")
    energies = []
    for line in file:
        if "Reaction coord              Energy" in line:
            for k in list(range( max_points(irc_file))):
                energy_k = (float(next(file).split()[2]))
                energies.append(energy_k)
    return energies



def coord_irc(irc_file):
    print('Analyzing each point of IRC ...')
    """Save the coordinate value of each points of IRC """
    file = open(irc_file,"r")
    coordinate = []
    for line in file:
        if "Reaction coord              Energy" in line:
            for k in list(range( max_points(irc_file) )):
                coord_k = (float(next(file).split()[1]))
                coordinate.append(coord_k)
        coor_zero_to_one = []
        for coord in coordinate:
            a_coord = round((coord - coordinate[0])/
                            (coordinate[max_points(irc_file) - 1] - coordinate[0]),4)
            coor_zero_to_one.append(a_coord)
    file.close()
    print('Done!')
    return coor_zero_to_one



def coord_transformation(irc_file, kr, kp, kts):
    """Perform a coordinate transformation from IRC to QRC """
    file = open(irc_file,"r")
    output_file = (open("sextic.dat", "w"))
    eact, erxn = [float(max(energy(irc_file))) ,
                  float(file.readlines()[-1].split()[2])]
    a, b, c, d, e = sexticrxn(eact, erxn, kr, kp, kts, "yes", output_file)
    print('Transforming coordinate from IRC to QRC ...')
    x = Symbol('x')
    init = Matrix([0.1])
    x_transformed = []
    for i in energy(irc_file):
        f = a * x ** 6 + b * x ** 5 + c * x ** 4 + d * x ** 3 + e * x ** 2 - i
        M = Matrix([f])
        JI = (M.jacobian([x]))**(-1)
        if i == 0:
            s0 = init
            iter = 0
            while(M.subs([(x,s0[0])]).norm()>(5e-8)):
                s0 = s0 - (M.subs([(x,s0[0])])*JI.subs([(x,s0[0])]))
                iter = iter + 1
                #print(s0[0])
            x_transformed.append(round(s0[0],4))
            sb = Matrix([s0[0] + 1/(max_points(irc_file) + 1)])
        elif i == erxn:
            while(M.subs([(x,sb[0])]).norm()>(5e-8)):
                sb = sb - (M.subs([(x,sb[0])])*JI.subs([(x,sb[0])]))
                iter = iter + 1
                #print(sb[0])
            x_transformed.append(round(sb[0],4))
            sb = Matrix([0.9999])
        else:
            while(M.subs([(x,sb[0])]).norm()>(5e-8)):
                sb = sb - (M.subs([(x,sb[0])])*JI.subs([(x,sb[0])]))
                iter = iter + 1
                #print(sb[0])
            x_transformed.append(round(sb[0],4))
            if sb[0] > 0.99:
                sb = Matrix([0.999])
            else:
                sb = Matrix([sb[0] + 1/(max_points(irc_file) + 1)])
    print('Done!')
    output_file.close()
    return x_transformed



def save_data(irc_file):
    """Save data of SRC analysis and the coordinate transformation """
    output_file = (open("sextic.dat", "a"))
    kr, kp = force_constant(irc_0, irc_1, irc_p_1, irc_p)
    kts = ts_force_constant(ts_file)
    irc = coord_irc(irc_file)
    src = coord_transformation(irc_file, kr, kp, kts)
    energy_of_each_point = energy(irc_file)
    sextic_plot(energy_of_each_point, irc, src, "plot", "save")
    output_file.write('\n=====================================================\n\n')
    output_file.write('\n              Coordinate transformation              \n')
    output_file.write('                   from IRC to SRC                   \n\n')
    output_file.write('\n=====================================================\n\n')
    output_file.write("\t IRC\t\t  SRC\t\t  Energy\n")
    output_file.write("\t\t\t  \t\t(kcal/mol)\n")
    output_file.write('=====================================================\n')
    print('Saving data ...')
    for point in list(range(max_points(irc_file))):
        output_file.write("\t" + str(irc[point]) + "\t\t" + str(src[point]) + "\t\t" +
                          str(energy_of_each_point[point])+ "\n")
    print('Done!')
    output_file.write('\n==========================END========================\n')
    output_file.close()


save_data(irc_file)
              

              
              
              
