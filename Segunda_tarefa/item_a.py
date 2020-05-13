#!/usr/bin/env python3
import numpy as np

def decompose_A(principal,sub):
    "Funcao que decompoe a matriz A em L e em D (A = LD(L^-1))"
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
    "Receive two square matrices and return a Matrix"
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
    "Receive a matrix and return it transposed"
    A = np.array(A) #make sure A is a np array not a list
    A_t = A.copy()

    for i in range(len(A)):
        for j in range(len(A[0])):
            A_t[j][i] = A[i][j]
    
    return A_t

def calculate_y(sub_L,b):
    "Receive sub diagonal of L(list) and the column matrix b and return y"
    # Ly = b

    y = np.array([b[0]])
    
    for i in range(len(b)-1):
        y = np.append(y, b[i+1] - sub_L[i]*y[i])
    
    return y

def calculate_z(diag_D,y):
    "Receive diagonal of D(list) and the column matrix y and return z"
    # Dz = y

    z = np.array([])
    
    for i in range(len(y)):
        z = np.append(z, y[i]/diag_D[i])
    
    return z

def calculate_x(sub_L,z):
    "Receive super diagonal of L transposed(list) and the column matrix z and return x"
    # L*(transposed)x = z"
    
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

def plot(us, letter, part, _u=None, erro=None):
    """
    Plot a graph using matplotlib
        us: array with heat_equation values (n=3)
        _u: array - y_utarget
        erro: list of floats
    """ 
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import sys
    import os

    # Save parameters
    current_path = os.path.abspath(__file__)
    current_path = current_path.split('/')
    current_path = current_path[:len(current_path) - 1]
    current_path = "/".join(current_path)

    # Figure parameters
    mpl.rcParams['lines.linewidth'] = 0.1
    plt.rcParams["figure.figsize"] = (20,2.2)
    
    fig, axs = plt.subplots(1,11, gridspec_kw={ 'hspace' : 1.5, 'wspace': 0.47}, constrained_layout=True)
    fig.suptitle('Plot para N = ' + str(len(us[0][0])-1))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    x_us = np.arange(0,1.0000000000001,1/(len(us[0][0])-1))
    us_dots = [8 for i in range(len(us[0][0]))] # list of dot sizes

    # Valores da solução exata
    x_utarget = np.arange(0,1,0.001)

    if(_u != None):
        y_target = np.array([_u(i) for i in x_utarget])
        
    target_dots = [0.2 for i in range(len(x_utarget))] # list of dot sizes
    
    # Plotting Graph
    for i in range(11):
        axs[i].scatter(x_us, us[0][i], s=us_dots, c='#119822')
        """ axs[0,i].set_xticks(np.arange(min(x_us), max(x_us)+1, 0.2)) """
        axs[i].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        axs.flat[i].yaxis.label.set_color('#119822')

    if(_u != None):
        axs.flat[10].scatter(x_utarget, y_target, s=target_dots, alpha=0.1)

    axs.flat[10].yaxis.set_label_position("right")
    axs.flat[10].yaxis.label.set_color('black')
    axs.flat[10].yaxis.label.set_fontsize(9)
    if(erro != None):
        axs.flat[10].set(ylabel="erro(T=1) = "+str(round(erro[0],10)))

    # save image as png
    if sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
        fig.savefig(r"Segunda_tarefa\figuras_{}\Figure of n = {}, parte {}.png".format(letter,len(us[0][0])-1, part), dpi=300)
    elif sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
        fig.savefig(current_path + "/figuras_{}" +"/Figure of n = {}, parte {}.png".format(letter,len(us[0][0])-1, part), dpi=300)
    else:
        print('--- AIX: saving fig at current directory ---')
        fig.savefig("letra_{}_figure of n = {}, parte {}.png".format(letter,len(us[0][0])-1, part), dpi=300)