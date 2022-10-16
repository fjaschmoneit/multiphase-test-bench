import numpy as np
from fieldAccess import *


def newField(shape):
    #shape = kwargs.get("shape")
    return np.ones(shape)



N = 5

shape = (N,N)


A = np.ones(shape)
for i in range(N):
    for j in range(N):
        A[i][j] = 10*i+j
print(A)




b = newField(shape)

b*=20


a = np.array(b)

a-=15

print(b)
print(a)

