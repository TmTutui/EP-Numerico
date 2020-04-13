#!/usr/bin/env python3
import numpy as np

def decompor_A(principal,sub):
    # Funcao que decompoe a matriz A em L e em D (A = LD(L^-1))
    # Recebe duas listas: a diagonal principal de A e a subdiagonal de A
    # Retorna duas listas: a diagonal Principal de D e a subdiagonal de L

    tamanho = len(principal)
    princ_D = np.array([])
    sub_L = np.array([])

    princ_D = np.append(princ_D,principal[0])


    for i in range(0,tamanho-1):
        l_i = sub[i]/princ_D[i] 

        sub_L = np.append(sub_L,l_i)

        princ_D = np.append(princ_D , principal[i+1]-(princ_D[i]*np.power(l_i,2)))

    return princ_D,sub_L

def dot_product(A,B): # work for two square matrices
    matrix = np.array([[0 for i in range(len(A))]])
    for i in range(len(A)-1):
        line = np.array([])
        for j in range(len(A)):
            line = np.append(line,0)

        matrix = np.append(matrix, [line], axis=0)


    for i in range(len(A)):  
        for j in range(len(B[0])):  
            for k in range(len(B)):  
                matrix[i][j] += A[i][k] * B[k][j] 

    return matrix

def transpose(A): # work for squace matrix
    A = np.array(A) #make sure A is a np array not a list
    A_t = A.copy()

    for i in range(len(A)):
        for j in range(len(A[0])):
            A_t[j][i] = A[i][j]
    
    return A_t



def main():


main()