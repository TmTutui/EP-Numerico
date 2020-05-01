import numpy as np
import math

def u(x):
    print(math.exp(1-x)*math.cos(5*x)) 
    print(np.exp(1-x)*np.cos(5*x))


""" print(max([0,1,2,3,4,5.5,5,6]))
print(np.max([0,1,2,3,4,5.5,5,6])) """

for i in range(11):
    u(i/10)

