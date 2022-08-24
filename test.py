import numpy as np


dim = 5

F = np.linspace(0,20,dim)


A = np.array([F,F], copy=False)


#A[:,0] = F
print(A)

F+=2
print(F)
print(A)