import numpy as np

a = np.zeros((2,6))
print(a)

pGrad = -2  #Pa/m
invCellDist = 10

a[:,:1] = 0+0.5*invCellDist*pGrad
for i in range(1,a.shape[1]):
    a[:,i:i+1] = a[:,i-1:i] + invCellDist*pGrad


print(a)