import numpy as np
from fieldAccess import *
import Interpolation
import Fields

nx = 4
ny = 1
a = Fields.newDataField((ny,nx))

for i in range(ny):
    for j in range(nx):
        a[i][j] = i+1 +10*(j+1)

print(a)

b = Interpolation.facesToNodes(a,'u')

print(b)
