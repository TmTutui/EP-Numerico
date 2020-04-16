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
    u_old = np.array([u0 for i in x_utarget])
    u_new = np.array([])
    
    for k in tqdm(range(1, M)):
        # adicionar u(k+1,0) na u_new
        u_new = np.array([g1])

        for i in range(1, N):
            u_new = np.append(u_new, u_old[i] + dt * ((u_old[i-1] - 2*u_old[i] + u_old[i+1]) / np.power(dx, 2) + _f(k*dt,i*dx) ))
        
        # adicionar u(k+1,N) na u_new
        u_new = np.append(u_new, g2)
        
        u_old = u_new.copy()
        
    # calcular o erro
    erro = np.max(abs(y_utarget-u_old))
        
    print('-'*15+'Heat Equation done'+'-'*15+'\n')
    return u_old, erro
    
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
    
    fig, axs = plt.subplots(3)
    fig.suptitle('Plot para N = ' + str(len(us[0])-1))
    plt.subplots_adjust(hspace = 0.4) # height space between subplots
    
    x_us = np.arange(0,1.0000000000001,1/(len(us[0])-1))
    us_dots = [2 for i in range(len(us[0]))] # list of dot sizes

    # Valores da solução exata
    x_utarget = np.arange(0,1,0.001)
    y_target = np.array([_u(x_utarget[i]) for i in range(len(x_utarget))])
    
    target_dots = [0.1 for i in range(len(y_target))] # list of dot sizes
    
    axs[0].scatter(x_us, us[0], s=us_dots)
    axs[0].scatter(x_utarget, y_target, s=target_dots, alpha=0.1)
    axs[0].set_title("u(t,x) | Lambda = 0.25 | erro(T=1) = "+str(erro[0]), 
                        fontdict={
                            'fontsize': 8,
                            'fontweight' : 0.3,
                        },
                        loc = 'center',
                    )

    axs[1].scatter(x_us, us[1], s=us_dots)
    axs[1].scatter(x_utarget, y_target, s=target_dots, alpha=0.1)
    axs[1].set_title("u(t,x) | Lambda = 0.50 | erro(T=1) = "+str(erro[1]), 
                        fontdict={
                            'fontsize': 8,
                            'fontweight' : 0.3,
                        },
                        loc = 'center',
                    )

    axs[2].scatter(x_us, us[2], s=us_dots)
    axs[2].scatter(x_utarget, y_target, s=target_dots, alpha=0.1)
    axs[2].set_title("u(t,x) | Lambda = 0.51 | erro(T=1) = "+str(erro[2]), 
                        fontdict={
                            'fontsize': 8,
                            'fontweight' : 0.3,
                        },
                        loc = 'center',
                    )

    # use same x label for every subplot
    for ax in fig.get_axes():
        ax.label_outer()

    # save image as png
    if sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
        fig.savefig(r"Primeira_tarefa\figuras_a\Figure of n = {}.png".format(len(us[0])-1), dpi=300)
    elif sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
        fig.savefig(current_path + "/figuras_a" +"/Figure of n = {}.png".format(len(us[0])-1), dpi=300)
    else:
        print('--- AIX: saving fig at current directory ---')
        fig.savefig("letra_a_figure of n = {}.png".format(len(us[0])-1), dpi=300)

def main():
    T = 1
    lamb_list = np.array([0.25 , 0.5 , 0.51])
    
    # u(0, x) = u0(x) em [0, 1]
    u0 = 0

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
        return 10*(np.power(x, 2))*(x - 1) - 60*x*t + 20*t 
    
    #solucao exata que precisamos nos aproximar:
    def _u(x):
        "Target solution"
        return 10*T*(np.power(x, 2))*(x - 1)        
    
    us = []
    erros = []
    
    for lamb in lamb_list:
        u_old, erro = heat_equation(u0, T, N, _f, lamb, g1, g2, _u)
        us.append(u_old)
        
        erros.append(erro)
    
    plot(us, _u, erros)
    print("--- %s seconds ---"%round(time.time() - start_time, 4))

main()