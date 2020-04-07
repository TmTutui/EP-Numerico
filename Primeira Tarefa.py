#!/usr/bin/env python3
import numpy as np    

def main():
    escolhido = False
    while (escolhido == False):
        alternative = input("Type an alternative (a,b ou c):")
        if(alternative.lower() == "a"):      
            escolhido = True
            letra_a()

        elif(alternative.lower() == "b"):
            escolhido = True
            letra_b()

        elif(alternative.lower() == "c"):
            escolhido = True
            letra_c()

        else:
            print(" You did not type an existing alternative! ")


def heat_equation(u0, T, _f, lamb):
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
    
    try:
        N = int(input("Type N: "))
    except:
        print("Wrong type! N must be integers!")
        N = int(input("Type N: "))
    
    delta_x = 1/N
    M = T*np.exp2(N)/lamb
    
    
    u_list = []
    """ 
    u_list = [
            [42,7,1287,325],
            [42,7,4,542],
            [42,7,981,325],
            [42,7,134,325]
            ]
    """

    for i in range(1, N, 1):
        u_list.append([])
        for k in range(0, M, 1):

            u_nextk_i = u_k_i + delta_t * ((u_k_lasti - 2*u_k_i + u_k_next_i) / np.exp2(delta_x) + _f(k*delta_t,i*delta_x))


def letra_a():
    T = 1
    lamb_list = [0.25 , 0.5 , 0.51]
    #u(0, x) = u0(x) em [0, 1]
    u0 = 0

    def _f(t, x):
        return 10*(np.exp2(x))*(x - 1) - 60*x*t + 20*t 
    
    #solucao exata que precisamos nos aproximar:
    def _u(t, x):
        return 10*t*(np.exp2(x))*(x - 1)

    for lamb in lamb_list:
        heat_equation(u0, T, _f, lamb)






