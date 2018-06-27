#!/usr/bin/python3

# -*- coding: utf-8 -*-

"""Sextic reaction coordinate
    
    This module contains the necessary functions to perform an optimization of the coefficients in the
    polinomial expansion of sextic reaction coordinate model (E = a*x**6 + b*x**5 + c*x**4 + d*x**3 + e*x**2).
    In order to perform this calculation, you must obtain the force constant of R, TS and P.
    To obtain them it is necessary to have a .fchk file from a frequency calculation performed in Gaussian
    for reactant and product, a .fchk file from a single point calculation of the second and the penultimate
    point of IRC and the .log of the ts optimization (with freq calculation).
    
    If you have problems with calculation, you can change the initial parameters. A list of previous tested
    initial parameters are showed below.
    
    """

from sympy import *
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.integrate import quad
from force_constant import *

DEFAULT_NAME = 'test'
WANT_PLT = 'no_save'

a, b, c, d, e = symbols('a, b, c, d, e')

#coef_matrix = Matrix([-640.011562944536, 2108.94580707296, -1788.09460758944, -202.381954261866, 527.432317722878]) # endo
#coef_matrix = Matrix([-698.0288158535812, 2000.4097266407093, -1657.308511379374, 30.403106250946294, 230.57449434129944]) # xts = 0.38
#coef_matrix = Matrix([-845.7973482679712, 2493.1951691064887, -1767.4200279106558, -561.5560584262729, 681.5782654984109]) # xts = 0.4972
#coef_matrix = Matrix([-588.785201131717, 1873.4092567113757, -1701.9579018291126, 152.84883805097772, 245.4750081984752]) # xts=0.4930
#coef_matrix = Matrix([-710.8704551  ,  2122.338051 ,   -1708.89548  ,  -126.1013714,    384.8759227]) # exo
#coef_matrix = Matrix([-675.0237837035049, 2201.994551010687, -2156.699474750905, 421.5304312838907, 201.18827615983392])


def initial_guess(erxn):
    """This function chooses the initial guess depending on the exothermicity
    of the reaction"""
    if erxn > 0:
        coef_matrix = Matrix([-640.0115629, 2108.945807, -1788.09460,
                              -202.38195426, 527.4323177])
    else:
        coef_matrix = Matrix([-698.0288158535812, 2000.4097266407093,
                              -1657.308511379374, 30.403106250946294,
                              230.57449434129944])
    return coef_matrix


def der(funcion,valor_a,valor_b,valor_c,valor_d,valor_e,dx):
    """This function calculates the numerical derivatives of each function"""
    f = Matrix([funcion])
    dfa = (f.subs([(a,(valor_a + dx)),(b,valor_b),(c,valor_c),(d,valor_d),(e,valor_e)]) -
           f.subs([(a,(valor_a - dx)),(b,valor_b),(c,valor_c),(d,valor_d),(e,valor_e)]))/(2*dx)
    dfb = (f.subs([(a,valor_a),(b,(valor_b + dx)),(c,valor_c),(d,valor_d),(e,valor_e)]) -
           f.subs([(a,valor_a),(b,(valor_b - dx)),(c,valor_c),(d,valor_d),(e,valor_e)]))/(2*dx)
    dfc = (f.subs([(a,valor_a),(b,valor_b),(c,(valor_c + dx)),(d,valor_d),(e,valor_e)]) -
           f.subs([(a,valor_a),(b,valor_b),(c,(valor_c - dx)),(d,valor_d),(e,valor_e)]))/(2*dx)
    dfd = (f.subs([(a,valor_a),(b,valor_b),(c,valor_c),(d,(valor_d + dx)),(e,valor_e)]) -
           f.subs([(a,valor_a),(b,valor_b),(c,valor_c),(d,(valor_d - dx)),(e,valor_e)]))/(2*dx)
    dfe = (f.subs([(a,valor_a),(b,valor_b),(c,valor_c),(d,valor_d),(e,(valor_e + dx))]) -
           f.subs([(a,valor_a),(b,valor_b),(c,valor_c),(d,valor_d),(e,(valor_e - dx))]))/(2*dx)
    return re(dfa[0].evalf()),re(dfb[0].evalf()),re(dfc[0].evalf()),re(dfd[0].evalf()),re(dfe[0].evalf())


def ts_search(coeffs):
    """This function perform a search of TS from the symbolic solution of a polynomium of sixth order.
        The result must complish: 0 < xts < 1 and to be a real number. It returns the TS position
        depending on the a, b, c, d and e coefficients of the expansion."""
    x = Symbol('x')
    f = a * x ** 6 + b * x ** 5 + c * x ** 4 + d * x ** 3 + e * x ** 2
    df = diff(f,x)
    sol = solve(df,x)
    solu = []
    for i in solve(df,x):
        if abs(im(i.subs([(a, coeffs[0]),(b, coeffs[1]),(c, coeffs[2]),
                          (d, coeffs[3]),(e, coeffs[4])]))) < (1e-10):
            solu.append(i)
    for i in solu:
        if (re(i.subs([(a, coeffs[0]),(b, coeffs[1]),(c, coeffs[2]),
                       (d, coeffs[3]),(e, coeffs[4])])).evalf() > 0.001 and
            re(i.subs([(a, coeffs[0]),(b, coeffs[1]),(c, coeffs[2]),
                       (d, coeffs[3]),(e, coeffs[4])])).evalf() < 0.999 and
            abs(im(i.subs([(a, coeffs[0]),(b, coeffs[1]),(c, coeffs[2]),
                           (d, coeffs[3]),(e, coeffs[4])]))).evalf() < 1e-10):
            xts = Matrix([i])
    return xts


def introduction(eact,erxn, kr, kp, kts, xts, coef_matrix, save_data, output_file):
    """ Generate the initial summary of data """
    if save_data=="yes":
        output_file.write('\n=====================================================\n\n\n')
        output_file.write('          SEXTIC REACTION COORDINATE ANALYSIS          \n\n')
        output_file.write('\n==================Initial Parameters=================\n\n\n')
        output_file.write('Eact         = '+  str(eact)+"\n")
        output_file.write('Erxn         = '+ str(erxn)+"\n")
        output_file.write('kr/kp        = '+ str(round((kr/kp),2))+"\n")
        output_file.write('kts/kr       = '+ str(round((kts/kr),2))+"\n")
        output_file.write('kts/kp       = '+ str(round((kts/kp),2))+"\n")
        output_file.write('ts_inicial   = '+ str(round(re(xts.subs([(a, coef_matrix[0]),(b, coef_matrix[1]),
                                           (c, coef_matrix[2]),(d, coef_matrix[3]),
                                           (e, coef_matrix[4])]))[0],4))+"\n")
    else:
        print('\n=====================================================\n\n')
        print('          SEXTIC REACTION COORDINATE ANALYSIS          \n')
        print('\n==================Initial Parameters=================\n')
        print('Eact         =', eact)
        print('Erxn         =', erxn)
        print('kr/kp        =', round((kr/kp),2))
        print('kts/kr       =', round((kts/kr),2))
        print('kts/kp       =', round((kts/kp),2))
        print('ts_inicial   = ',round(re(xts.subs([(a, coef_matrix[0]),(b, coef_matrix[1]),
                                                   (c, coef_matrix[2]),(d, coef_matrix[3]),
                                                   (e, coef_matrix[4])]))[0],4))


def crit_convergence(functions, coef_matrix):
    """ Evaluate the convergence criteria """
    convergence_crit = functions.subs([(a, coef_matrix[0]),(b, coef_matrix[1]),
                                       (c, coef_matrix[2]),(d, coef_matrix[3]),
                                       (e, coef_matrix[4])]).norm()
    return convergence_crit


def jacob(functions, coef_matrix):
    """ Obtain the jacobian matrix """
    jacobiano = Matrix([der(functions[0],coef_matrix[0], coef_matrix[1],
                            coef_matrix[2], coef_matrix[3], coef_matrix[4], 1e-5),
                        der(functions[1],coef_matrix[0], coef_matrix[1],
                            coef_matrix[2],coef_matrix[3], coef_matrix[4], 1e-5),
                        der(functions[2],coef_matrix[0], coef_matrix[1],
                            coef_matrix[2],coef_matrix[3], coef_matrix[4], 1e-5),
                        der(functions[3],coef_matrix[0], coef_matrix[1],
                            coef_matrix[2],coef_matrix[3], coef_matrix[4], 1e-5),
                        der(functions[4],coef_matrix[0], coef_matrix[1],
                            coef_matrix[2],coef_matrix[3], coef_matrix[4], 1e-5)])
    return jacobiano


def final_results(opt_coeff, xts, convergence, save_data, output_file, j):
    """ Generate the final summary of data """
    ts_position = (re(xts.subs([(a, opt_coeff[0]),(b, opt_coeff[1]),
                                (c, opt_coeff[2]),(d, opt_coeff[3]),
                                (e, opt_coeff[4])]))[0])
    kr,kp,kts = force_constant(opt_coeff, ts_position)
    arc_len = arc_lenght(opt_coeff)
    cp1, cp2, w1, w2, w3, w4 = reaction_force(opt_coeff, ts_position)
    if save_data=="yes":
        output_file.write('\n====================Final Results====================\n\n')
        output_file.write('Convergence reached after '+str((j+1))+' cycles\n')
        output_file.write('Optimized coef.      = '+str(opt_coeff)+'\n')
        output_file.write('Final Convergence    = '+ str(convergence)+'\n')
        output_file.write('Final TS position    = '+str(round(ts_position,2))+'\n')
        output_file.write('\n=====================================================\n\n')
        output_file.write('Force_constant_r     = '+str(round(kr,2))+'\n')
        output_file.write('Force_constant_ts    = '+str(round(kts,2))+'\n')
        output_file.write('Force_constant_p     = '+str(round(kp,2))+'\n')
        output_file.write('\n=====================================================\n\n')
        output_file.write('Arc lenght           = '+ str(round(arc_len[0],2))+'\n')
        output_file.write('\n===============Reaction Force Analysis===============\n\n')
        output_file.write('F_min                = '+ str(round(cp1,4))+'\n')
        output_file.write('F_max                = '+ str(round(cp2,4))+'\n\n')
        output_file.write('W1                   = '+ str(round(w1,2))+ ' ('+str(round(w1*100/(w1+w2),2))+'% )'+'\n')
        output_file.write('W2                   = '+ str(round(w2,2))+ ' ('+str(round(w2*100/(w1+w2),2))+'% )'+'\n')
        output_file.write('W3                   = '+ str(round(w3,2))+ ' ('+str(round(w3*100/(w3+w4),2))+'% )'+'\n')
        output_file.write('W4                   = '+ str(round(w4,2))+ ' ('+str(round(w4*100/(w3+w4),2))+'% )'+'\n')
        output_file.write('\n=====================================================\n')
    else:
        print('\n====================Final Results====================\n')
        print('Convergence reached after',(j+1),'cycles')
        print('Optimized coef.      =',opt_coeff)
        print('Final Convergence    =', convergence)
        print('Final TS position =',round(ts_position,2))
        print('\n=====================================================\n')
        print('Force_constant_r     =',round(kr,2))
        print('Force_constant_ts    =',round(kts,2))
        print('Force_constant_p     =',round(kp,2))
        print('\n=====================================================\n')
        print('Arc lenght           =', round(arc_len[0],2))
        print('\n===============Reaction Force Analysis===============\n')
        print('F_min                =', round(cp1,4))
        print('F_max                =', round(cp2,4),'\n')
        print('W1                   =', round(w1,2), '(',round(w1*100/(w1+w2),2),'% )')
        print('W2                   =', round(w2,2), '(',round(w2*100/(w1+w2),2),'% )')
        print('W3                   =', round(w3,2), '(',round(w3*100/(w3+w4),2),'% )')
        print('W4                   =', round(w4,2), '(',round(w4*100/(w3+w4),2),'% )')
        print('\n=====================================================\n')


def force_constant(opt_coeff, xts):
    """ Evaluate the second derivative in R, TS and P """
    kr = 2 * opt_coeff[4]
    kp = (30 * opt_coeff[0] + 20 * opt_coeff[1] + 12 * opt_coeff[2] + 6 *
          opt_coeff[3] + 2 * opt_coeff[4])
    kts = (30 * opt_coeff[0] * xts ** 4 + 20 * opt_coeff[1] * xts ** 3 + 12 *
           opt_coeff[2] * xts ** 2 + 6 * opt_coeff[3] * xts + 2 * opt_coeff[4])
    return kr, kp, kts


def reaction_force(opt_coef, ts_position):
    """ Perform a reaction force analysis """
    x = Symbol('x')
    y = (opt_coef[0] * x ** 6 + opt_coef[1] * x ** 5 + opt_coef[2] * x ** 4
                      + opt_coef[3] * x ** 3 + opt_coef[4] * x ** 2)
    dy = diff(y,x)
    cp_f = solve(diff(dy,x))
    cp_real=[]
    for i in cp_f:
        if abs(im(i)) < (1e-5) and re(i) > 0 and re(i) < 1:
            cp_real.append(abs(i))
    cp_real.sort()
    w1 = integrate(dy,(x,0,cp_real[0]))
    w2 = integrate(dy,(x,cp_real[0],ts_position))
    w3 = integrate(dy,(x,ts_position, cp_real[1]))
    w4 = integrate(dy,(x,cp_real[1],1))
    return cp_real[0], cp_real[1], w1, w2, w3, w4


def arc_lenght(coef):
    """ Calculate the arc lenght of the sextic function between 0 and 1 """
    f = lambda x:(sqrt((6 * coef[0] * x ** 5 + 5 * coef[1] * x ** 4 +
                        4 * coef[2] * x ** 3+ 3 * coef[3] * x ** 2 +
                        2 * coef[4] * x) ** 2 + 1))
    return quad(f,0,1)


# -------Sextic Reaction Coordinate-------


def sexticrxn(eact, erxn, kr, kp, kts, save_data, output_file):
    """This function performs an iterative search of the optimized coefficients
        of the sextic expansion, using the newton's method to solve a non-linear equation"""
    coef_matrix = initial_guess(erxn)
    xts = ts_search(coef_matrix)
    introduction(eact, erxn, kr, kp, kts, xts, coef_matrix, save_data, output_file)
    x = Symbol('x')
    func1 = a + b + c + d + e - erxn
    func2 = (a * re(xts[0]) ** 6 + b * re(xts[0]) ** 5 + c * re(xts[0]) ** 4
             + d * re(xts[0]) ** 3 + e * re(xts[0]) ** 2 - eact)
    func3 = 6 * a + 5 * b + 4 * c + 3 * d + 2 * e
    func4 = ((30 * a * re(xts[0]) ** 4 + 20 * b * re(xts[0]) ** 3 + 12 * c *
              re(xts[0]) ** 2 + 6 * d *re(xts[0]) + 2 * e) * kr - (2 * e) * kts)
    func5 = ((30 * a + 20 * b  + 12 * c  + 6 * d + 2 * e) * kts -
             (30 * a * re(xts[0]) ** 4 + 20 * b * re(xts[0]) ** 3 + 12 * c *
              re(xts[0]) ** 2 + 6 * d *re(xts[0]) + 2 * e) * kp)
    print('Constructing Function\'s matrix ...')
    functions_matrix = Matrix([func1, func2, func3, func4, func5])
    print('Done!')
    print('Constructing Inverted Jacobian matrix ...')
    jac_inv = jacob(functions_matrix, coef_matrix) ** (-1)
    print('Done!')
    print('\n=====================================================\n')
    j = 0
    convergence = crit_convergence(functions_matrix, coef_matrix)
    print('\nIterative process begins ...\n')
    while (convergence > 1e-8):
        print('Cycle',(j+1))
        ts_position = (re(xts.subs([(a, coef_matrix[0]),(b, coef_matrix[1]),
                                    (c, coef_matrix[2]),(d, coef_matrix[3]),
                                    (e, coef_matrix[4])]))[0])
        print('xts=', ts_position)
        print('Convergence criteria =', convergence)
        coef_matrix = coef_matrix - (jac_inv *
                                     functions_matrix.subs([(a, coef_matrix[0]),(b, coef_matrix[1]),
                                                            (c, coef_matrix[2]),(d, coef_matrix[3]),
                                                            (e, coef_matrix[4])]))
        xts = ts_search(coef_matrix)
        print('Updating and Inverting Jacobian ...')
        jac_inv = jacob(functions_matrix, coef_matrix) ** (-1)
        print('Done!')
        print(coef_matrix,'\n')
        convergence = crit_convergence(functions_matrix, coef_matrix)
        j = j + 1
    opt_coef = [float(coef_matrix[0]), float(coef_matrix[1]),
                float(coef_matrix[2]),float(coef_matrix[3]),
                float(coef_matrix[4])]
    final_results(opt_coef, xts, convergence, save_data, output_file, j)
    return opt_coef


def sextic_plot(e_irc, coord_irc, coord_src, name=DEFAULT_NAME,want_plot=WANT_PLT):
    rcParams['font.family'] = 'Times New Roman'
    plt.plot(coord_irc, e_irc, label='IRC', marker='o', fillstyle='none',
                 markersize='2.5', linewidth=0.5, c='g')
    plt.plot(coord_src, e_irc, label='SRC', marker='o', fillstyle='none',
             markersize='3', linewidth=0.5, c='orange')
    plt.title('IRC vs. SRC', fontsize=13)
    plt.ylabel('Potential Energy', fontsize=13)
    plt.xlabel('Reaction coordinate', fontsize=13)
    plt.legend()
    if want_plot == 'save':
        return plt.savefig(name+'.png', dpi=1080)
    else:
        plt.show()
    return






