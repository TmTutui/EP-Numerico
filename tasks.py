#!/usr/bin/env python3
import os

import numpy as np
from tqdm import tqdm

from Segunda_tarefa.item_a import plot, decompose_A, calculate_x, calculate_y, calculate_z

def Task01(N, pontos, k):
    " Calcula uk(T,x)"
    nf = len(pontos)
    
    def _u0(x):
        "Initial distribution"
        return 0
    def _g1(t):
        "Distribution on x = 0"
        return 0
    def _g2(t):
        "Distribution on x =N"
        return 0
    def r(t):
        return 10*(1 + np.cos(5*t))

    def _ghk(x):
        h = 1/N 
        pk = pontos[k-1]
        if(pk - h/2 <= x <= pk + h/2):
            return 1/h
        else: 
            return 0
    
    def _f(t, x, nf):
        return r(t) * _ghk(x)
            
    def heat_equation(_u0, _f, _g1, _g2, T=1, N = N):
        """
        Heat Equation:
            u0: uo(x) - math function
            N: int (input)
            T: float
            i = 1, ..., N-1
            k = 0, ..., M-1
            f: math function - f(t,x)
            u: Heat Equation - u(t, x)
            xi = i∆x, i = 0, · · · , N, com ∆x = 1/N. Para a discretização temporal definimos ∆t = T /M, e
            calculamos aproximações nos instantes tk = k∆t, k = 1, · · · , M. 
            A variável u(t, x) descreve a temperatura no instante t na posição x, sendo a distribuição inicial u0(x) dada

        return: 
            u_new = u(T,x)
        """
        
        print('-'*15+'Heat Equation in progress'+'-'*15+'\n')
        
        dx = 1/N
        dt = dx
        M = int(T/dt)
        lamb = 1/dx

        # used in u exata
        x_utarget = np.arange(0, 1.0000000001, dx)

        u_old = np.array([_u0(i) for i in x_utarget])

        # matrix A
        A_diag = np.array([(1+lamb) for i in range(N-1)])
        A_sub = np.array([(-lamb/2) for i in range(N-2)])

        diag_D, sub_L = decompose_A(A_diag,A_sub)

        # Ax = b ou seja A*u_new[1:N-1] = b
        for k in tqdm(range(0, M)):
            # adicionar u(k+1,0) na u_new
            u_new = np.array([_g1((k+1)*dt)])

            # create b 
            b = np.array([])
            b = np.append(b, u_old[1] + (lamb/2)*(_g1((k+1)*dt) + u_old[0] - 2*u_old[1] + u_old[2]) + (dt/2)*(_f(dt*(k+1),dx*1) + _f(dt*k,dx*1)))
            for i in range(2, N-1):
                b = np.append(b, u_old[i] + (lamb/2)*(u_old[i-1] - 2*u_old[i] + u_old[i+1]) + (dt/2)*(_f(dt*(k+1),dx*i) + _f(dt*k,dx*i)))
            b = np.append(b, u_old[N-1] + (lamb/2)*(_g2((k+1)*dt) + u_old[N-2] - 2*u_old[N-1] + u_old[N]) + (dt/2)*(_f(dt*(k+1),dx*(N-1)) + _f(dt*k,dx*(N-1))))

            # find x
            y = calculate_y(sub_L,b)
            z = calculate_z(diag_D,y)
            x = calculate_x(sub_L,z)

            for x_element in x:
                u_new = np.append(u_new, x_element)
            
            # adicionar u(k+1,N) na u_new
            u_new = np.append(u_new, _g2((k+1)*dt))
            
            u_old = u_new.copy()
            # print(u_old)

        return u_new
    
    uk_T = heat_equation(_u0, _f, _g1, _g2)
    return uk_T[1,len(uk_T)-2]

def Task02(uT, uks):

    """Item b do exercício
        uT: 
        us: multi-dimensional array of array
    """
    
    def _inner_product(vetor1, vetor2):
        """
            vetor1: numpy array
            vetor2: numpy array
        """
        return sum(vetor1[:]*vetor2[:])
    
    def _matriz(uks):
        """
        uks: 2D numpy array 
        """
        A = np.array([]) # 2D matrix with dot product of uks
        b = np.array([]) # right hand side 1D matrix
        
        for k1 in range(len(uks)):
            A = np.append(A, [_inner_product(uks[k2], uks[k1]) for k2 in range(len(uks))])
            b = np.append(b, _inner_product(uT, uks[k1]))

        # [ ] TODO: ver se não tem que adicionar vetor x ao cálculo do sistem normal
        return A, b
        
        for k1 in range(len(uks)):
            for k2 in range(len(uks)):
                A = np.append(A, _inner_product(uks[k1], uks[k2]))

def Task03(M_sim):
    def find_D_and_L(M_sim):
        D = np.array([])
        
        m_dimension = len(M_sim)
        L = np.array([np.zeros(m_dimension)])
        L[0,0] = 1
        
        for i in range(1,m_dimension):
            L = np.append(L, [np.zeros(m_dimension)], axis=0)
            L[i,i] = 1

        for line in range(m_dimension):
            for column in range(m_dimension):
                if(line >= column):
                    if(line==column):
                        sum = 0
                        for n in range(line):
                            sum+=D[n]*np.square(L[line][n])
                        D = np.append(D, M_sim[line][column] - sum)
                    
                    else:
                        sum = 0
                        for n in range(column):
                            sum+= L[line][n]*L[column][n]*D[n]
                        L[line,column] = (M_sim[line][column] - sum)/D[column]
        
        return D, L

    def Teste():
        m_tst = [
            [4,12,-16],
            [12,37,-43],
            [-16,-43,98]
        ]
    
        D,L = find_D_and_L(m_tst)

        print(D)
        print(L)

    Teste()
                


Task03([0,0,0])