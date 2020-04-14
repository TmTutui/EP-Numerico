#!/usr/bin/env python3
import numpy as np

def decompose_A(principal,sub):
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

@DeprecationWarning
def dot_product(A,B): # work for two square matrices
    # Receive two square matrices and return a Matrix
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

@DeprecationWarning
def transpose(A): # work for squace matrix
    # Receive a matrix and return it transposed
    A = np.array(A) #make sure A is a np array not a list
    A_t = A.copy()

    for i in range(len(A)):
        for j in range(len(A[0])):
            A_t[j][i] = A[i][j]
    
    return A_t

def calculate_y(sub_L,b):
    # Ly = b
    # Receive sub diagonal of L(list) and the column matrix b and return y
    y = np.array([b[0]])
    
    for i in range(len(b)-1):
        y = np.append(y, b[i+1] - sub_L[i]*y[i])
    
    return y

def calculate_z(diag_D,y):
    # Dz = y
    # Receive diagonal of D(list) and the column matrix y and return z
    z = np.array([])
    
    for i in range(len(y)):
        z = np.append(z, y[i]/diag_D[i])
    
    return z

def calculate_x(sub_L,z):
    # L*(transposed)x = z
    # Receive super diagonal of L transposed(list) and the column matrix z and return x
    x = np.zeros(len(z))
    x[len(x)-1] = z[len(z)-1]
    
    for i in range(len(z)-2,-1,-1):
        x[i] = z[i] - sub_L[i]*x[i+1]
    
    return x

def teste4():
    A = [
        [32.5,65,0,0],
        [65,131,5,0],
        [0,5,28,6],
        [0,0,6,16],
    ]
    b=[
       5677.75,
       11477.3,
       23791.08,
       61804.16
    ]

    diag_D, sub_L = decompose_A([32.5,131,28,16],[65,5,6])
    # Ly = b, Dz = y e L*x = z
    y = calculate_y(sub_L,b)
    

    z = calculate_z(diag_D,y)
    x = calculate_x(sub_L,z)

    print(x)


