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

def newCellField(mesh, value):
    cf = np.ndarray(shape=(mesh.cells_y, mesh.cells_x), dtype=float, order='C')
    setInitialValue(cf, value)
    return cf

def newFaceField_x(mesh, value):
    ff = np.ndarray(shape=(mesh.cells_y, mesh.cells_x+1), dtype=float, order='C')
    setInitialValue(ff, value)
    return ff

def newFaceField_y(mesh, value):
    ff = np.ndarray(shape=(mesh.cells_y+1, mesh.cells_x), dtype=float, order='C')
    setInitialValue(ff, value)
    return ff

    #
# class newFaceField_xold(object):
#     def __init__(self, mesh, value=0):
#         self.mesh = mesh
#         self.shape = (mesh.cells_y, mesh.cells_x+1)
#         self.nbFaces = self.shape[0]*self.shape[1]
#         self.field = np.ndarray(shape=self.shape, dtype=float, order='C')
#         self.internal = self.field[:,1:-1]
#         self.boundary = self.field[:,::self.mesh.cells_x]
#         self.setInitialValue(value)

    # def setInternalValues(self, field):
    #     self.field[:,1:-1] = field

    def setBoundaryValues_w(self, field):
        self.field[:,:1] = field

    def setBoundaryValues_e(self, field):
        self.field[:, -1:] = field

    # def setInitialValue(self, value):
    #     self.field.fill(value)
    #
    # def fillWithConsecutiveValues(self):
    #     i = 0
    #     for y in range(len(self.field)):
    #         for x in range(len(self.field[0])):
    #             self.field[y][x] = i
    #             i += 1

# class faceField_y(object):
#     def __init__(self, mesh, value=0):
#         self.mesh = mesh
#
#         self.shape = (mesh.cells_y+1, mesh.cells_x)
#         self.nbFaces = self.shape[0] * self.shape[1]
#         self.field = np.ndarray(shape=self.shape, dtype=float, order='C')
#         self.internal = self.field[1:-1, :]
#         self.boundary = self.field[::mesh.cells_y, :]
#         self.setInitialValue(value)
#
#     def setInternalValues(self, field):
#         self.field[1:-1,:] = field
#
#     def setInitialValue(self, value):
#         self.field.fill(value)

    def setBoundaryValues_n(self, field):
        self.field[0:1, :] = field

    def setBoundaryValues_s(self, field):
        self.field[-1:, :] = field

    # def fillWithConsecutiveValues(self):
    #     i = 0
    #     for y in range(len(self.field)):
    #         for x in range(len(self.field[0])):
    #             self.field[y][x] = i
    #             i += 1


