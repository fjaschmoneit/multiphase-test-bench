# primitive Fields hold data

import numpy as np

def setInitialValue(primitiveField, value):
    primitiveField.fill(value)

def fillWithConsecutiveValues(primitiveField):
    i = 0
    for y in range(len(primitiveField)):
        for x in range(len(primitiveField[0])):
            primitiveField[y][x] = i
            i += 1

def fillWithRandomIntegers(primitiveField):
    primitiveField[:,:] = np.random.randint( low=0, high=10, size=primitiveField.shape )

#needed??
def newVertexField(mesh, value):
    vf = np.ndarray(shape=(mesh._cells_y+1, mesh._cells_x+1), dtype=float, order='C'  )
    setInitialValue(vf, value)
    return vf

def newCellField(mesh, value):
    cf = np.ndarray(shape=(mesh._cells_y, mesh._cells_x), dtype=float, order='C')
    setInitialValue(cf, value)
    return cf

# neededs??
def newFaceField_x(mesh, value):
    ff = np.ndarray(shape=(mesh._cells_y, mesh._cells_x+1), dtype=float, order='C')
    setInitialValue(ff, value)
    return ff

def newFaceField_y(mesh, value):
    ff = np.ndarray(shape=(mesh._cells_y+1, mesh._cells_x), dtype=float, order='C')
    setInitialValue(ff, value)
    return ff
