import time
start_time = time.time()
import numpy as np
from tqdm import tqdm
import os
import sys

from item_a import decompose_A, calculate_x,calculate_y,calculate_z

current_path = os.path.abspath(__file__)
current_path = current_path.split('/')
current_path = current_path[:len(current_path) - 1]
current_path = "/".join(current_path)

def heat_equation(_u0, T, N, _f, lamb, _g1, _g2, _u):
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
    dt = dx
    M = int(T/dt)     

    # used in u exata
    x_utarget = np.arange(0, 1.0000000001, dx)
    y_utarget = np.array([_u(x_utarget[i]) for i in range(len(x_utarget))])

    # used in aprox
    u_old = np.array([_u0(i) for i in x_utarget])

    # matrix A
    A_diag = np.array([(1+2*lamb) for i in range(N-1)])
    A_sub = np.array([(-lamb) for i in range(N-2)])

    diag_D, sub_L = decompose_A(A_diag,A_sub)

    # Ax = b ou seja A*u_new[1:N-1] = b
    for k in tqdm(range(0, M)):
        # adicionar u(k+1,0) na u_new
        u_new = np.array([_g1(0)])

        # create b 
        b = np.array([])
        for i in range(1, N):
            # it is possible to do everything in a loop cause g1=g2=0
            b = np.append(b, u_old[i] + dt*_f(dt*(k+1),dx*i))

        # find x
        y = calculate_y(sub_L,b)
        z = calculate_z(diag_D,y)
        x = calculate_x(sub_L,z)

        for x_element in x:
            u_new = np.append(u_new, x_element)
        
        # adicionar u(k+1,N) na u_new
        u_new = np.append(u_new, _g2(0))
        """ print(u_new) """
        
        u_old = u_new.copy()
        print(u_old)

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
        fig.savefig(r"Segunda_tarefa\figuras_b\Figure of n = {}.png".format(len(us[0])-1), dpi=300)
    elif sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
        fig.savefig(current_path + "/figuras_b" +"/Figure of n = {}.png".format(len(us[0])-1), dpi=300)
    else:
        print('--- AIX: saving fig at current directory ---')
        fig.savefig("letra_b_figure of n = {}.png".format(len(us[0])-1), dpi=300)

def main():
    T = 1
    lamb_list = np.array([0.25, 0.5, 0.51])
    
    def _f(t, x):
        "Descrição da fonte de calor ao longo do tempo"
        return 10*np.cos(10*t) * x**2 * (1-x)**2 - (1 + np.sin(10*t))*(12*x**2 - 12*x + 2)
    
    def _u0(x):
        "Condição de contorno"
        return np.power(x, 2) * np.power((1 - x), 2)

    #solucao exata que precisamos nos aproximar:
    def _u(x):
        "Target solution"
        return (1 + np.sin(10*1)) * x**2 * (1 - x)**2   
    
    def _g1(t):
        "Condição de fronteira x = 0."
        return 0
    
    def _g2(t):
        "Condição de fronteira x = 1."
        return 0
    
    try:
        N = int(input("Type N: "))
    except:
        print("Wrong type! N must be an integer!")
        N = int(input("Type N: "))
    
    us = []
    erros = []

    for lamb in lamb_list:
        u_old, erro = heat_equation(_u0, T, N, _f, lamb, _g1, _g2, _u)
        us.append(u_old)
        
        erros.append(erro)
    
    plot(us, _u, erros)
    print("--- %s seconds ---"%round(time.time() - start_time, 4))

main()
