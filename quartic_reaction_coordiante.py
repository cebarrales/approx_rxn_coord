#!/usr/bin/python3

# -*- coding: utf-8 -*-

"""Quartic reaction coordinate
    
    This module contains a function, which allows to perform an optimization of the coefficients
    in the polinomial expansion of quartic reaction coordinate model (E = a*x**4 + b*x**3 + c*x**2).
    Also, it contains the main function to perfomr an additional analysis once the optimal coefficients
    are obtained. You just need to give the activation energy and reaction energy in the quarticrxn
    function.
    
    Example:
    
    a,b,c = quarticrxn(20.0, -12.0)
    quartic_plot1(a,b,c, 'yes')
    
    It will be perform a calculation of the optimized coefficients showing the results in the screen.
    It also will plot the data, saving a file called "test_plot.png". If you want to change the name,
    you must specify as the fifth parameter of the quartic_plot function without the extension png.
    
    """

from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.integrate import quad


DEFAULT_NAME = 'test_plot'
WANT_PLT = 'no'
SAVE = 'no'
OUT = 'test_data'
a, b, c = symbols('a, b, c')


def der(funcion,value_a,value_b,value_c,dx):
    """ Calculate the numerical derivative of each function """
    f = Matrix([funcion])
    dfa = (f.subs([(a,(value_a + dx)),(b,value_b),(c,value_c)]) -
           f.subs([(a,(value_a - dx)),(b,value_b),(c,value_c)]))/(2*dx)
    dfb = (f.subs([(a,value_a),(b,(value_b + dx)),(c,value_c)]) -
           f.subs([(a,value_a),(b,(value_b - dx)),(c,value_c)]))/(2*dx)
    dfc = (f.subs([(a,value_a),(b,value_b),(c,(value_c + dx))]) -
           f.subs([(a,value_a),(b,value_b),(c,(value_c - dx))]))/(2*dx)
    return dfa[0],dfb[0],dfc[0]


def jacob(functions, coef_matrix):
    """ Obtain the jacobian matrix """
    jacobiano = Matrix([der(functions[0],coef_matrix[0], coef_matrix[1],
                            coef_matrix[2], 1e-5),
                        der(functions[1],coef_matrix[0], coef_matrix[1],
                            coef_matrix[2], 1e-5),
                        der(functions[2],coef_matrix[0], coef_matrix[1],
                            coef_matrix[2], 1e-5),
                        ])
    return jacobiano


def crit_convergence(functions, coef_matrix):
    """ Evaluate the convergence criteria """
    convergence_crit = functions.subs([(a, coef_matrix[0]),(b, coef_matrix[1]),
                                       (c, coef_matrix[2])]).norm()
    return convergence_crit


def introduction(eact,erxn, save_data, output_file):
    """ Generate the initial summary of data """
    if save_data=="yes":
        output_file.write('=====================================================\n\n\n')
        output_file.write('       QUARTIC REACTION COORDINATE ANALYSIS      \n\n')
        output_file.write('\n==================Initial Parameters=================\n\n\n')
        output_file.write("Eact \t= " + str(eact) + "\n")
        output_file.write("Erxn \t= " + str(erxn) + "\n")
    else:
        print('=====================================================\n\n')
        print('       QUARTIC REACTION COORDINATE ANALYSIS      \n')
        print('\n==================Initial Parameters=================\n')
        print('Eact =', eact)
        print('Erxn =', erxn)


def final_results(opt_coeff, ts_position, convergence, save_data, output_file, j):
    """ Generate the final summary of data """
    kr,kp,kts = force_constant(opt_coeff, ts_position)
    arc_len = arc_lenght(opt_coeff)
    cp1, cp2, w1, w2, w3, w4 = reaction_force(opt_coeff, ts_position)
    if save_data=="yes":
        output_file.write('\n====================Final Results====================\n\n')
        output_file.write("Convergence reached after "+str((j+1))+" cycles\n")
        output_file.write("Optimized coef.      = "+str(opt_coeff)+"\n")
        output_file.write("Final Convergence    = "+ str(convergence)+"\n")
        output_file.write("Final TS position    = "+str(round(ts_position,2))+"\n")
        output_file.write("\n=====================================================\n\n")
        output_file.write("Force_constant_r     = "+str(round(kr,2))+"\n")
        output_file.write("Force_constant_ts    = "+str(round(kts,2))+"\n")
        output_file.write("Force_constant_p     = "+str(round(kp,2))+"\n")
        output_file.write('\n=====================================================\n\n')
        output_file.write("Arc lenght           = "+ str(round(arc_len[0],2))+"\n")
        output_file.write('\n===============Reaction Force Analysis===============\n\n')
        output_file.write("F_min                = "+ str(round(cp1,4))+"\n")
        output_file.write("F_max                = "+ str(round(cp2,4))+"\n")
        output_file.write("Width of TS_region   = "+ str(round(cp2-cp1,4))+"\n\n")
        output_file.write("W1                   = "+ str(round(w1,2))+ "("+str(round(w1*100/(w1+w2),2))+"% )\n")
        output_file.write("W2                   = "+ str(round(w2,2))+ "("+str(round(w2*100/(w1+w2),2))+"% )\n")
        output_file.write("W3                   = "+ str(round(w3,2))+ "("+str(round(w3*100/(w3+w4),2))+"% )\n")
        output_file.write("W4                   = "+ str(round(w4,2))+ "("+str(round(w4*100/(w3+w4),2))+"% )\n")
        output_file.write('\n=====================================================')
    else:
        print('\n====================Final Results====================\n')
        print('Convergence reached after',(j+1),'cycles')
        print('Optimized coef.      =',opt_coeff)
        print('Final Convergence    =', convergence)
        print('Final TS position    =',round(ts_position,2))
        print('\n=====================================================\n')
        print('Force_constant_r     =',round(kr,2))
        print('Force_constant_ts    =',round(kts,2))
        print('Force_constant_p     =',round(kp,2))
        print('\n=====================================================\n')
        print('Arc lenght           =', round(arc_len[0],2))
        print('\n===============Reaction Force Analysis===============\n')
        print('F_min                =', round(cp1,4))
        print('F_max                =', round(cp2,4))
        print('Width of TS_region   =', round(cp2-cp1,4),'\n')
        print('W1                   =', round(w1,2), '(',round(w1*100/(w1+w2),2),'% )')
        print('W2                   =', round(w2,2), '(',round(w2*100/(w1+w2),2),'% )')
        print('W3                   =', round(w3,2), '(',round(w3*100/(w3+w4),2),'% )')
        print('W4                   =', round(w4,2), '(',round(w4*100/(w3+w4),2),'% )')
        print('\n=====================================================\n')


def reaction_force(opt_coeff, xts):
    """ Perform a reaction force analysis """
    x = Symbol('x')
    y = opt_coeff[0] * x ** 4 + opt_coeff[1] * x ** 3 + opt_coeff[2] * x ** 2
    cp_f = solve(diff(diff(y)))
    work1 = integrate(diff(y), (x, 0, cp_f[0]))
    work2 = integrate(diff(y), (x, cp_f[0], xts))
    work3 = integrate(diff(y), (x, xts, cp_f[1]))
    work4 = integrate(diff(y), (x, cp_f[1], 1))
    return cp_f[0], cp_f[1], work1, work2, work3, work4


def force_constant(opt_coeff, xts):
    """ Evaluate the second derivative in R, TS and P """
    kr = 2 * opt_coeff[2]
    kp = 12 * opt_coeff[0] + 6 * opt_coeff[1] + 2 * opt_coeff[2]
    kts = 12 * opt_coeff[0] * xts ** 2 + 6 * opt_coeff[1] * xts + 2 * opt_coeff[2]
    return kr, kp, kts


def arc_lenght(coef):
    """ Calculate the arc lenght of the quartic function between 0 and 1 """
    f = lambda x:(sqrt((4 * coef[0] * x ** 3+ 3 * coef[1] * x ** 2 + 2 * coef[2] * x) ** 2 + 1))
    return quad(f,0,1)


# -------Quartic Reaction Coordinate-------
# ----------Iterative process--------------

def quarticrxn(eact, erxn, save_data=SAVE, output_file=OUT):
    """Perform a coefficient optimization of the quartic reaction coordinate.
        
        Parameters
        ----------
        eact: float
        Activation energy of the chemical reaction (E_transition_state - E_reactant).
        erxn: float
        Reaction energy of the chemical reaction (E_product - E_reactant).
        
        Returns
        -------
        The optimized functions will be in a list (a, b, c).
        
        """
    introduction(eact, erxn, save_data, output_file)
    x = Symbol('x')
    xts = ((-6 * b) / (8 * a) - 1)
    coef_matrix = Matrix([100, -200, 200])
    func1 = a + b + c - erxn
    func2 = 4 * a + 3 * b + 2 * c
    func3 = a * xts ** 4 + b * xts ** 3 + c * xts ** 2 - eact
    functions_matrix = Matrix([func1, func2, func3])
    jacobiano = Matrix([der(func1,coef_matrix[0], coef_matrix[1], coef_matrix[2], 1e-5),
                        der(func2,coef_matrix[0], coef_matrix[1], coef_matrix[2], 1e-5),
                        der(func3,coef_matrix[0], coef_matrix[1], coef_matrix[2], 1e-5)])
    jac_inv = jacobiano ** -1
    j = 0
    convergence = crit_convergence(functions_matrix, coef_matrix)
    while (convergence > 1e-8):
        coef_matrix = coef_matrix - (jac_inv * functions_matrix.subs([(a, coef_matrix[0]),
                                                                      (b, coef_matrix[1]),
                                                                      (c, coef_matrix[2])]))
        jac_inv = jacob(functions_matrix, coef_matrix) ** (-1)
        convergence = crit_convergence(functions_matrix, coef_matrix)
        j = j + 1
    opt_coeff = [float(coef_matrix[0]), float(coef_matrix[1]), float(coef_matrix[2])]
    ts_position = ((-6 * opt_coeff[1]) / (8 * opt_coeff[0]) - 1)
    kr = 2 * opt_coeff[2]
    kp = 12 * opt_coeff[0] + 6 * opt_coeff[1] + 2 * opt_coeff[2]
    final_results(opt_coeff, ts_position, convergence, save_data, output_file, j)
    return opt_coeff


# ----------Plot Final Results--------------


def quartic_plot1(a, b, c, want_plot=WANT_PLT, name=DEFAULT_NAME):
    """ Collect the data and plot it """
    def frange(ini, fin, delta):
        j = ini
        while j < fin:
            yield j
            j += delta
    xaxis = list(frange(0, 1.01, 0.01))
    yaxis = []
    for n in frange(0, 1.01, 0.01):
        yaxis.append(a * n ** 4
                     + b * n ** 3
                     + c * n ** 2
                     )
    rcParams['font.family'] = 'Times New Roman'
    plt.plot(xaxis, yaxis, marker='o', fillstyle='none',
             markersize='3', linewidth=0.5)
    plt.ylabel('Potential Energy', fontsize=13)
    plt.xlabel('Quartic Reaction coordinate', fontsize=13)
    if want_plot == 'yes':
        return plt.savefig(name+'.png', dpi=1080)
    else:
        return plt.show()



def quartic_plot(e_irc, coord_irc, coord_qrc, name=DEFAULT_NAME,want_plot=WANT_PLT):
    """ Collect the data from coordinate transformation and plot the QRC and IRC """
    rcParams['font.family'] = 'Times New Roman'
    plt.plot(coord_irc, e_irc, label='IRC', marker='o', fillstyle='none',
             markersize='2.5', linewidth=0.5, c='g')
    plt.plot(coord_qrc, e_irc, label='QRC', marker='o', fillstyle='none',
                      markersize='3', linewidth=0.5)
    plt.title('IRC vs. QRC', fontsize=13)
    plt.ylabel('Potential Energy (kcal/mol)', fontsize=13)
    plt.xlabel('Reaction coordinate', fontsize=13)
    plt.legend(fontsize=13)
    if want_plot == 'save':
        return plt.savefig(name+'.png', dpi=1080)
    else:
        plt.show()
    return


#a,b,c = quarticrxn(20.0, -12.0)
#quartic_plot1(a,b,c, 'yes')

