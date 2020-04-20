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

def heat_equation(u0, T, N, _f, lamb, g1, g2, _u):
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
    y_utarget = np.array([_u(x_utarget[i]) for i in range(len(x_utarget))])


    # used in aprox
    u_old = np.array([u0(i) for i in x_utarget])

    # u for every 0.1 units of time
    u_interval = np.array([u_old])
    list_times = [i for i in range(0, M +1 ,M//10)]
    
    for k in tqdm(range(0, M)):
        # adicionar u(k+1,0) na u_new
        u_new = np.array([g1(k*dt)])

        for i in range(1, N):
            u_new = np.append(u_new, u_old[i] + dt * ((u_old[i-1] - 2*u_old[i] + u_old[i+1]) / np.power(dx, 2) + _f(k*dt,i*dx)))
        
        # adicionar u(k+1,N) na u_new
        u_new = np.append(u_new, g2(k*dt))
        
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

    Save figures at figuras_b
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
    y_target = np.array([_u(x_utarget[i]) for i in range(len(x_utarget))])
    
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
        fig.savefig(r"Primeira_tarefa\figuras_b\Figure of n = {}.png".format(len(us[0][0])-1), dpi=300)
    elif sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
        fig.savefig(current_path + "/figuras_b" +"/Figure of n = {}.png".format(len(us[0][0])-1), dpi=300)
    else:
        print('--- AIX: saving fig at current directory ---')
        fig.savefig("letra_b_figure of n = {}.png".format(len(us[0][0])-1), dpi=300)

def main():
    
    def _u0(x):
        "Distribuição inicial."
        return np.exp(-x)
    
    def _g1(t):
        "Condição de fronteira x = 0."
        return np.exp(t)
    
    def _g2(t):
        "Condição de fronteira x = 1."
        return np.exp(t-1)*np.cos(5*t)
    
    def _f(t, x):
        "Descrição da fonte de calor ao longo do tempo"
        return np.exp(t-x)*5*(5*np.power(t,2)*np.cos(5*t*x) - (x + 2*t)*np.sin(5*t*x))
        
    T = 1
    
    try:
        N = int(input("Type N: "))
    except:
        print("Wrong type! N must be an integer!")
        N = int(input("Type N: "))
        
    def _u(x):
        "Target solution (with t=T=1)"
        return np.exp(1-x)*np.cos(5*x)    

    us = []
    erros = []
    
    for lamb in np.array([0.25 , 0.5 , 0.51]):
        u_old, erro = heat_equation(_u0, T, N, _f, lamb, _g1, _g2, _u)
        us.append(u_old)
        
        erros.append(erro)
        
    print('Erro: \n\n', erros, '\n\n')    
    plot(us, _u, erros)
    print("--- %s seconds ---"%round(time.time() - start_time, 4))

main()