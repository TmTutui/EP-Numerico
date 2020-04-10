#!/usr/bin/env python3
import time
start_time = time.time()
import numpy as np   
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

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
    """
    
    print('-'*15+'Heat Equation in progress'+'-'*15+'\n')
    
    dx = 1/N
    M = int(T*np.power(N, 2)/lamb)
    dt = T/M    

    # used in aprox
    u_old = np.arange(0, 1.0000000001, dx)
    u_new = np.array([])
    
    for k in tqdm(range(1, M)):
        # adicionar u(k+1,0) na u_new
        u_new = np.append(u_new, 0)

        for i in range(1, N):
            u_new = np.append(u_new, u_old[i] + dt * ((u_old[i-1] - 2*u_old[i] + u_old[i+1]) / np.power(dx, 2) + _f(k*dt,i*dx,dx)))
        
        # adicionar u(k+1,N) na u_new
        u_new = np.append(u_new, 0)
        
        u_old = u_new.copy()
        u_new = []
        
    print('-'*15+'Heat Equation done'+'-'*15+'\n')
    return u_old


def plot(us):
    """
    Plot a graph using matplotlib
        us: array with heat_equation values (n=3)

    Save figures at figuras_c
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['lines.linewidth'] = 0.1
    
    fig, axs = plt.subplots(3,11)
    fig.suptitle('Plot para N = ' + str(len(us[0])-1))
    plt.subplots_adjust(hspace = 0.4) # height space between subplots
    
    x_us = np.arange(0,1.0000000000001,1/(len(us[0])-1))
    us_dots = [2 for i in range(len(us[0]))] # list of dot sizes
    
    
    axs[0,0].scatter(x_us, us[0], s=us_dots)
    axs.flat[0].set(xlabel='tempo = ' + str(0), ylabel='Lambda = 0.25')


    axs[1,0].scatter(x_us, us[1], s=us_dots)
    axs.flat[1].set(xlabel='tempo = ' + str(0.1), ylabel='Lambda = 0.5')


    axs[2,0].scatter(x_us, us[2], s=us_dots)
    axs.flat[2].set(xlabel='tempo = ' + str(0.2), ylabel='Lambda = 0.51')

    # use same x label for every subplot
    for ax in fig.get_axes():
        ax.label_outer()

    # save image as png
    try:
        fig.savefig(r"Primeira_tarefa\figuras_c\Figure of n = {}.png".format(len(us[0])-1), dpi=300)
    except:
        fig.savefig(r"Primeira_tarefa/figuras_c/Figure of n = {}.png".format(len(us[0])-1), dpi=300)

def main():
    
    def _u0(x):
        "Distribuição inicial."
        return 0
    
    def _g1(t):
        "Condição de fronteira x = 0."
        return 0
    
    def _g2(t):
        "Condição de fronteira x = 1."
        return 0
    
    def _f(t, x, dx):
        "Descrição da fonte de calor ao longo do tempo, = r(t) * Gh(x)  "
        p = 0.25
        # h = dx

        if (p-dx <= x <= p):
            "gh(x) poderia assumir o valor 1/h em p e variar linearmente de 0 a 1/h no intervalo [p − h, p]"
            return 10000*(1-2*np.power(t,2)) * ((1/np.power(dx,2))*(x + dx - p))

        if (p < x <= p + dx):
            "e (gh(x) poderia assumir o valor) de 1/h a 0 no intervalo [p, p + h], sendo nula no restante do dom´ınio."
            return 10000*(1-2*np.power(t,2)) * ((1/np.power(dx,2))*(-x + dx + p))
        
        return 0
        
                
    T = 1
    lamb_list = [0.25 , 0.5 , 0.51]
    
    try:
        N = int(input("Type N: "))
    except:
        print("Wrong type! N must be an integer!")
        N = int(input("Type N: "))
        
    def _u(x):
        "Target solution."
        return np.exp(T-x)*np.cos(5*T*x)    

    us = []
    
    for lamb in lamb_list:
        u_old = heat_equation(_u0, T, N, _f, lamb, _g1, _g2, _u)
        us.append(u_old)

        
    plot(us)
    print("--- %s seconds ---"%round(time.time() - start_time, 4))

main()