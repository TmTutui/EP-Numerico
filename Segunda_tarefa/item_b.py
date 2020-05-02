import time
start_time = time.time()
import numpy as np
from tqdm import tqdm

from Segunda_tarefa.item_a import plot, decompose_A, calculate_x, calculate_y, calculate_z

def main():
    try:
        N = int(input("Type N: "))
    except:
        print("Wrong type! N must be an integer!")
        N = int(input("Type N: "))

    part_a(N)
    part_b(N)
    part_c(N)

def heat_equation(_u0, T, N, _f, _g1, _g2, _u=None):
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
    lamb = 1/dx

    # used in u exata
    x_utarget = np.arange(0, 1.0000000001, dx)
    if(_u != None):
        y_utarget = np.array([_u(i) for i in x_utarget])

    # used in aprox
    u_old = np.array([_u0(i) for i in x_utarget])

    # u for every 0.1 units of time
    u_interval = np.array([u_old])
    list_times = [i for i in range(0, M +1 ,M//10)]

    # matrix A
    A_diag = np.array([(1+2*lamb) for i in range(N-1)])
    A_sub = np.array([(-lamb) for i in range(N-2)])

    diag_D, sub_L = decompose_A(A_diag,A_sub)

    # Ax = b ou seja A*u_new[1:N-1] = b
    for k in tqdm(range(0, M)):
        # adicionar u(k+1,0) na u_new
        u_new = np.array([_g1((k+1)*dt)])

        # create b 
        b = np.array([u_old[1] + dt*_f(dt*(k+1),dx) + lamb*_g1(dt*(k+1))])
        for i in range(2, N-1):
            # it is possible to do everything in a loop cause g1=g2=0
            b = np.append(b, u_old[i] + dt*_f(dt*(k+1),dx*i))

        b = np.append(b, u_old[N-1] + dt*_f(dt*(k+1),dx*(N-1)) + lamb*_g2(dt*(k+1)) )

        # find x
        y = calculate_y(sub_L,b)
        z = calculate_z(diag_D,y)
        x = calculate_x(sub_L,z)

        for x_element in x:
            u_new = np.append(u_new, x_element)
        
        # adicionar u(k+1,N) na u_new
        u_new = np.append(u_new, _g2((k+1)*dt))
        """ print(u_new) """
        
        u_old = u_new.copy()

        if( (k+1) in list_times ):
            u_interval = np.append(u_interval, [u_old], axis = 0)

    if(_u != None):
        # calcular o erro
        erro = np.max(abs(y_utarget-u_old))
            
        print('-'*15+'Heat Equation done'+'-'*15+'\n')
        return u_interval, erro

    else:
        print('-'*15+'Heat Equation done'+'-'*15+'\n')
        return u_interval

def part_a(N):
    T = 1
    
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
    
    us = []
    erros = []

    
    u_old, erro = heat_equation(_u0, T, N, _f, _g1, _g2, _u)
    us.append(u_old)
    
    erros.append(erro)
    
    plot(us, "b", "A", _u, erros)

def part_b(N):
    
    def _u0(x):
        "Distribuição inicial."
        return np.exp(-x)
    
    def _g1(t):
        "Condição de fronteira x = 0."
        return np.exp(t)
    
    def _g2(t):
        "Condição de fronteira x = 1."
        return np.exp(t-1)*np.cos(5*t)
    
    def _f(t, x):
        "Descrição da fonte de calor ao longo do tempo"
        return np.exp(t-x)*5*(5*np.power(t,2)*np.cos(5*t*x) - (x + 2*t)*np.sin(5*t*x))
        # return -np.exp(t-x)*(5*x*np.sin(5*t*x) - np.cos(5*t*x) + 10*t*np.sin(5*t*x) + (1-25*t*t)*np.cos(5*t*x) )
        
    T = 1
        
    def _u(x):
        "Target solution (with t=T=1)"
        return np.exp(1-x)*np.cos(5*x)    

    us = []
    erros = []
    
    
    u_old, erro = heat_equation(_u0, T, N, _f, _g1, _g2, _u)
    us.append(u_old)
    
    erros.append(erro)
    
    plot(us, 'b', "B", _u, erros)


def part_c(N):
    
    def _u0(x):
        "Distribuição inicial."
        return 0
    
    def _g1(t):
        "Condição de fronteira x = 0."
        return 0
    
    def _g2(t):
        "Condição de fronteira x = 1."
        return 0
    
    def _f(t, x):
        "Descrição da fonte de calor ao longo do tempo, = r(t) * Gh(x)  "
        p = 0.25
        # h = dx
        dx = 1/N

        if (p-dx <= x <= p):
            "gh(x) poderia assumir o valor 1/h em p e variar linearmente de 0 a 1/h no intervalo [p − h, p]"
            return 10000*(1-2*np.power(t,2)) * ((1/np.power(dx,2))*(x + dx - p))

        elif (p < x <= p + dx):
            "e (gh(x) poderia assumir o valor de 1/h a 0 no intervalo [p, p + h], sendo nula no restante do domínio."
            return 10000*(1-2*np.power(t,2)) * ((1/np.power(dx,2))*(-x + dx + p))
        
        else:
            return 0
                 
    T = 1
    us = []
    
    u_olds = heat_equation(_u0, T, N, _f, _g1, _g2)
    us.append(u_olds)

        
    plot(us, 'b', "C")

main()