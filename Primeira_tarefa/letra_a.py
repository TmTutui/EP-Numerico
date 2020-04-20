#!/usr/bin/env python3
import time
start_time = time.time()
import sys
import numpy as np
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")
import os

current_path = os.path.abspath(__file__)
current_path = current_path.split('/')
current_path = current_path[:len(current_path) - 1]
current_path = "/".join(current_path)
    
def heat_equation(_u0, T, N, _f, lamb, g1, g2, _u):
    """
    Heat Equation:
        u0: uo(x) - math function
        N: int (input)
        M: int (input)
        T: float
        i = 1, ..., N-1
        k = 0, ..., M-1
        f: math function - f(t,x)
        u: Heat Equation - u(t, x)
        xi = i∆x, i = 0, · · · , N, com ∆x = 1/N. Para a discretização temporal definimos ∆t = T /M, e
        calculamos aproximações nos instantes tk = k∆t, k = 1, · · · , M. 
        A variável u(t, x) descreve a temperatura no instante t na posição x, sendo a distribuição inicial u0(x) dada

    return: 
        u_old: array
        erro: list
    """
    
    print('-'*15+'Heat Equation in progress'+'-'*15+'\n')
    
    dx = 1/N
    M = int(T*np.power(N, 2)/lamb)
    dt = T/M    

    # used in u exata
    x_utarget = np.arange(0, 1.0000000001, dx)
    y_utarget = np.array([_u(T, x_utarget[i]) for i in range(len(x_utarget))])

    # used in aprox
    u_old = np.array([_u0(x_utarget[i]) for i in range(len(x_utarget))])

    # u for every 0.1 units of time
    u_interval = np.array([u_old])
    list_times = [i for i in range(0, M +1 ,M//10)]
    
    for k in tqdm(range(0, M)):
        # adicionar u(k+1,0) na u_new
        u_new = np.array([g1])

        for i in range(1, N):
            u_new = np.append(u_new, u_old[i] + dt * ((u_old[i-1] - 2*u_old[i] + u_old[i+1]) / np.power(dx, 2) + _f(k*dt,i*dx) ))   

        # adicionar u(k+1,N) na u_new
        u_new = np.append(u_new, g2)

        u_old = u_new.copy()

        if( (k+1) in list_times ):
            u_interval = np.append(u_interval, [u_old], axis = 0)
    
    # calcular o erro
    erro = np.max(abs(y_utarget-u_old))
        
    print('-'*15+'Heat Equation done'+'-'*15+'\n')
    return u_interval, erro   
    
def plot(us, _u, erro):
    """
    Plot a graph using matplotlib
        us: array with heat_equation values (n=3)
        _u: array - y_utarget
        erro: list of floats

    Save figures at figuras_a
    """ 
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['lines.linewidth'] = 0.1
    plt.rcParams["figure.figsize"] = (40,5)
    
    fig, axs = plt.subplots(3,11, gridspec_kw={ 'hspace' : 0.45, 'wspace': 0.47})
    fig.suptitle('Plot para N = ' + str(len(us[0][0])-1))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    x_us = np.arange(0,1.0000000000001,1/(len(us[0][0])-1))
    us_dots = [2 for i in range(len(us[0][0]))] # list of dot sizes

    # Valores da solução exata
    x_utarget = np.arange(0,1,0.001)
    y_target = np.array([_u(1, x_utarget[i]) for i in range(len(x_utarget))])
    
    target_dots = [0.1 for i in range(len(y_target))] # list of dot sizes
    

    for i in range(11):
        axs[0,i].scatter(x_us, us[0][i], s=us_dots, c='#80ab4e')
        """ axs[0,i].set_xticks(np.arange(min(x_us), max(x_us)+1, 0.2)) """
        axs[0,i].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axs.flat[i].yaxis.label.set_color('#80ab4e')

        axs[1,i].scatter(x_us, us[1][i], s=us_dots, c='#FF8C00')
        """ axs[1,i].set_xticks(np.arange(min(x_us), max(x_us)+1, 0.2)) """
        axs[1,i].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axs.flat[i + 11].yaxis.label.set_color('#FF8C00')

        axs[2,i].scatter(x_us, us[2][i], s=us_dots, c='red')
        """ axs[2,i].set_xticks(np.arange(min(x_us), max(x_us)+1, 0.2)) """
        axs[2,i].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axs.flat[i + 22].set(xlabel='tempo = ' + str(i/10))
        axs.flat[i + 22].yaxis.label.set_color('red')


    axs.flat[0].set(ylabel='Lambda = 0.25')
    axs.flat[11].set(ylabel='Lambda = 0.5')
    axs.flat[22].set(ylabel='Lambda = 0.51')


    axs.flat[10].scatter(x_utarget, y_target, s=target_dots, alpha=0.1)
    axs.flat[21].scatter(x_utarget, y_target, s=target_dots, alpha=0.1)
    axs.flat[32].scatter(x_utarget, y_target, s=target_dots, alpha=0.1)

    axs.flat[10].yaxis.set_label_position("right")
    axs.flat[10].yaxis.label.set_color('black')
    axs.flat[10].yaxis.label.set_fontsize(6)
    axs.flat[10].set(ylabel="erro(T=1) = "+str(round(erro[0],10)))
    
    axs.flat[21].yaxis.set_label_position("right")
    axs.flat[21].yaxis.label.set_color('black')
    axs.flat[21].yaxis.label.set_fontsize(6)
    axs.flat[21].set(ylabel="erro(T=1) = "+str(round(erro[1],10)))

    axs.flat[32].yaxis.set_label_position("right")
    axs.flat[32].yaxis.label.set_color('black')
    axs.flat[32].yaxis.label.set_fontsize(6)
    axs.flat[32].set(ylabel="erro(T=1) = "+str(round(erro[2],10)))


    # save image as png
    if sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
        fig.savefig(r"Primeira_tarefa\figuras_a\Figure of n = {}.png".format(len(us[0][0])-1), dpi=300)
    elif sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
        fig.savefig(current_path + "/figuras_a" +"/Figure of n = {}.png".format(len(us[0][0])-1), dpi=300)
    else:
        print('--- AIX: saving fig at current directory ---')
        fig.savefig("letra_a_figure of n = {}.png".format(len(us[0][0])-1), dpi=300)

def main():
    T = 1
    
    # condições de fronteira nulas
    g1 = 0
    g2 = 0
    
    try:
        N = int(input("Type N: "))
    except:
        print("Wrong type! N must be an integer!")
        N = int(input("Type N: "))

    def _f(t, x):
        "Descrição da fonte de calor ao longo do tempo"
        return 10*np.cos(10*t) * x**2 * (1-x)**2 - (1 + np.sin(10*t))*(12*x**2 - 12*x + 2)
    
    def _u0(x):
        "Condição de contorno"
        return np.power(x, 2) * np.power((1 - x), 2)

    #solucao exata que precisamos nos aproximar:
    def _u(t, x):
        "Target solution"
        return (1 + np.sin(10*t)) * x**2 * (1 - x)**2     
    
    us = []
    erros = []
    
    for lamb in np.array([0.25 , 0.5 , 0.51]):
        u_olds, erro = heat_equation(_u0, T, N, _f, lamb, g1, g2, _u)
        us.append(u_olds)
        
        erros.append(erro)
    
    plot(us, _u, erros)
    print("--- %s seconds ---"%round(time.time() - start_time, 4))

main()