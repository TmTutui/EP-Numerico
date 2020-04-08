#!/usr/bin/env python3
import numpy as np   
import math 
from tqdm import tqdm

def main():
    escolhido = False
    
    while (escolhido == False):
        alternative = input("Type an alternative (a,b ou c):")
        if(alternative.lower() == "a"):      
            escolhido = True
            letra_a()

        # elif(alternative.lower() == "b"):
        #     escolhido = True
        #     letra_b()

        # elif(alternative.lower() == "c"):
        #     escolhido = True
        #     letra_c()

        else:
            print(" You did not type an existing alternative! ")


def heat_equation(u0, T, N, _f, lamb, g1, g2):
    """
    Fórmula de diferenças finitas:
        N: int (input)
        M: int (input)
        T: float
        i = 1, ..., N-1
        k = 0, ..., M-1
        f: math function - f(t,x)
        u: Heat Equation - u(t, x)


        xi = i∆x, i = 0, · · · , N, com ∆x = 1/N. Para a discretiza¸c˜ao temporal definimos ∆t = T /M, e
        calculamos aproxima¸c˜oes nos instantes tk = k∆t, k = 1, · · · , M. Denotamos a aproxima¸c˜ao para a
        solu¸c˜ao nos pontos de malha u(tk, xi) por u_k_i.
        
        A variável u(t, x) descreve a temperatura no instante t na posi¸c˜ao x, sendo a distribuição inicial u0(x) dada
    """
    
    print('-'*15+'Heat Equation in progress'+'-'*15+'\n')
    
    dx = 1/N
    M = int(T*np.power(N, 2)/lamb)
    dt = T/M    
    
    u_old = np.array([u0 for i in range(N+1)])
    type(u_old)
    
    u_new = np.array([])

    
    for k in tqdm(range(1, M)):
        # adicionar u(k+1,0) na u_new
        u_new = np.append(u_new, g1)

        for i in range(1, N):
            u_new = np.append(u_new, u_old[i] + dt * ((u_old[i-1] - 2*u_old[i] + u_old[i+1]) / np.power(dx, 2) + _f(k*dt,i*dx)))
        
        # adicionar u(k+1,N) na u_new
        u_new = np.append(u_new, g2)

        u_old = u_new.copy()
        u_new = []        

    
    print('-'*15+'Heat Equation done'+'-'*15+'\n')
    return u_old
    
def plot(us, _u):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['lines.linewidth'] = 0.1
    
    fig = plt.figure()
    
    ax_lamb1 = fig.add_subplot(4, 4, 1)
    ax_lamb2 = fig.add_subplot(4, 4, 2)
    ax_lamb3 = fig.add_subplot(4, 4, 3)

    print(us[0])
    print(us[1])
    print(us[2]) 
    
    def dot_size_us(n_dots):
        return [1 for i in range(n_dots)]


    x_us = np.arange(0,1.0000000000001,1/(len(us[0])-1))
    us_dots = [2 for i in range(len(us[0]))] #list of dot sizes

    #Valores da solucao exata
    x_utarget = np.arange(0,1,0.001)
    y_target = np.array([_u(x_utarget[i]) for i in range(len(x_utarget))])
    target_dots = [0.1 for i in range(len(y_target))] #list of dot sizes

    
    ax_lamb1.scatter(x_us, us[0], s=us_dots)
    ax_lamb1.scatter(x_utarget, y_target, s=target_dots, alpha=0.1)

    ax_lamb2.scatter(x_us, us[1], s=us_dots)
    ax_lamb2.scatter(x_utarget, y_target, s=target_dots, alpha=0.1)

    ax_lamb3.scatter(x_us, us[2], s=us_dots)
    ax_lamb3.scatter(x_utarget, y_target, s=target_dots, alpha=0.1)

    #save image as png
    fig.savefig("figura.png", dpi=300)

def letra_a():
    T = 1
    lamb_list = [0.25 , 0.5 , 0.51]
    #u(0, x) = u0(x) em [0, 1]
    u0 = 0

    #condicoes de fronteira nulas
    g1 = 0
    g2 = 0

    try:
        N = int(input("Type N: "))
    except:
        print("Wrong type! N must be an integer!")
        N = int(input("Type N: "))


    def _f(t, x):
        return 10*(np.power(x, 2))*(x - 1) - 60*x*t + 20*t 
    
    #solucao exata que precisamos nos aproximar:
    def _u(x):
        return 10*T*(np.power(x, 2))*(x - 1)

    us = []
    
    for lamb in lamb_list:
        u_old = heat_equation(u0, T, N, _f, lamb, g1, g2)
        us.append(u_old)
        
    plot(us, _u)


main()