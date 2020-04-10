#!/usr/bin/env python3
import time
start_time = time.time()
import numpy as np   
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

def heat_equation(T, N, _f, lamb):
    """
    Heat Equation:
        u0: int = 0
        N: int (input)
        M: int (input)
        T: float
        i = 1, ..., N-1
        k = 0, ..., M-1
        f: math function - f(t,x)
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

    # used in aprox, u0 = 0
    u_old = np.array([0 for i in range(N+1)])
    u_new = np.array(np.array([]))

    # u for every 0.1 units of time
    u_interval = np.array([u_old])
    list_times = [i for i in range(0, M +1 ,M//10)]
    
    for k in tqdm(range(1, M)):
        # adicionar u(k+1,0) = 0 na u_new
        u_new = np.append(u_new, 0)

        for i in range(1, N):
            u_new = np.append(u_new, u_old[i] + dt * ((u_old[i-1] - 2*u_old[i] + u_old[i+1]) / np.power(dx, 2) + _f(k*dt,i*dx,dx)))
        
        # adicionar u(k+1,N) = 0 na u_new
        u_new = np.append(u_new, 0)
        
        u_old = u_new.copy()
        u_new = []

        if( (k+1) in list_times ):
            u_interval = np.append(u_interval, [u_old], axis = 0)

    print(u_interval)
        
    print('-'*15+'Heat Equation done'+'-'*15+'\n')
    return u_interval


def plot(us):
    """
    Plot a graph using matplotlib
        us: array with heat_equation values (n=3)

    Save figures at figuras_c
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

    # save image as png
    try:
        fig.savefig(r"Primeira_tarefa\figuras_c\Figure of n = {}.png".format(len(us[0][0])-1), dpi=300)
    except:
        fig.savefig(r"Primeira_tarefa/figuras_c/Figure of n = {}.png".format(len(us[0][0])-1), dpi=300)

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
        u_olds = heat_equation(T, N, _f, lamb)
        us.append(u_olds)

        
    plot(us)
    print("--- %s seconds ---"%round(time.time() - start_time, 4))

main()