import importlib
import numpy as np

import ScalarField
importlib.reload(ScalarField)


class hectorField:
    # access functions to underlying primitive structure
#    def __init__(self, mesh, value=0, primitiveFieldX=None, primitiveFieldY=None):
    def __init__(self, u, v):
        print("moin")

        # these are scalar fields:
        self._u = u
        self._v = v

        # these are ndarrays
        self._e = self._u.e
        self._w = self._u.w
        self._n = self._v.n
        self._s = self._v.s

        self._be = self._u.be
        self._bw = self._u.bw
        self._bn = self._v.bn
        self._bs = self._v.bs

    @classmethod
    def fromMesh(cls, mesh, value=0.0):
        return cls( *cls.UVfromMesh( mesh, value) )

    @staticmethod
    def UVfromMesh(mesh, value = 0):
        u = ScalarField.scalarField.fromShape(mesh.shapeStaggered_U, value)
        v = ScalarField.scalarField.fromShape(mesh.shapeStaggered_V, value)
        return (u, v)


#-------- defining algebra
    def __add__(self, other):
        u = self._u + other._u
        v = self._v + other._v
        return vectorField(u,v)

    def __mul__(self, other):
        if isinstance(other, vectorField):
            return vectorField(self._u * other._u, self._v * other._v)
        elif isinstance(other, type(1.0)):
            return vectorField(self._u * other, self._v * other)
        else:
            print("multiplication not defined")

    #-------- defining all getters and setters:

    @property
    def u(self):
        return self._u
    @u.setter
    def u(self, x):
        self._u.data[:, :] = x
    @property
    def v(self):
        return self._v
    @v.setter
    def v(self, x):
        self._v.data[:, :] = x

    # @property
    # def internalEntries_NS(self):
    #     return self._v_internal
    # @internalEntries_NS.setter
    # def internalEntries_NS(self, x):
    #     self._v_internal[:, :] = x
    #
    # @property
    # def internalEntries_EW(self):
    #     return self._u_internal
    # @internalEntries_EW.setter
    # def internalEntries_EW(self, x):
    #     self._u_internal[:, :] = x
    #
    @property
    def e(self):
        return self._e
    @e.setter
    def e(self, x):
        self._e[:,:] = x

    @property
    def w(self):
        return self._w
    @w.setter
    def w(self, x):
        self._w[:,:] = x

    @property
    def n(self):
        return self._n
    @n.setter
    def n(self, x):
        self._n[:,:] = x

    @property
    def s(self):
        return self._s
    @s.setter
    def s(self, x):
        self._s[:,:] = x

    @property
    def be(self):
        return self._be
    @be.setter
    def be(self, x):
        self._be = x

    @property
    def bw(self):
        return self._bw
    @bw.setter
    def bw(self, x):
        self._bw = x

    @property
    def bn(self):
        return self._bn
    @bn.setter
    def bn(self, x):
        self._bn = x

    @property
    def bs(self):
        return self._bs
    @bs.setter
    def bs(self, x):
        self._bs = x











class staggeredFluxField_U:

    def __init__(self, mesh, value=0):
        self._mesh = mesh

        self._u = np.ndarray(shape=(mesh._cells_y, mesh._cells_x), dtype=float, order='C')
        self._v = np.ndarray(shape=(mesh._cells_y+1, mesh._cells_x-1), dtype=float, order='C')

        self._u_internal = self._u[:, 1:-1]
        self._v_internal = self._v[1:-1, :]

        self._e = self._u[:, 1:]
        self._w = self._u[:, :-1]
        self._n = self._v[:-1, :]
        self._s = self._v[1:, :]

        self._be = self._u[:, -1:]
        self._bw = self._u[:, :1]
        self._bn = self._v[:1, :]
        self._bs = self._v[-1:, :]

    @property
    def entries_NS(self):
        return self._v

    @entries_NS.setter
    def entries_NS(self, x):
        self._v[:, :] = x

    @property
    def entries_EW(self):
        return self._u

    @entries_EW.setter
    def entries_EW(self, x):
        self._u[:, :] = x

    @property
    def internalEntries_NS(self):
        return self._v_internal

    @internalEntries_NS.setter
    def internalEntries_NS(self, x):
        self._v_internal[:, :] = x

    @property
    def internalEntries_EW(self):
        return self._u_internal

    @internalEntries_EW.setter
    def internalEntries_EW(self, x):
        self._u_internal[:, :] = x

    @property
    def e(self):
        return self._e

    @e.setter
    def e(self, x):
        self._e[:, :] = x

    @property
    def w(self):
        return self._w

    @w.setter
    def w(self, x):
        self._w[:, :] = x

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, x):
        self._n[:, :] = x

    @property
    def s(self):
        return self._s

    @s.setter
    def s(self, x):
        self._s[:, :] = x

    @property
    def be(self):
        return self._be

    @be.setter
    def be(self, x):
        self._be[:, :] = x

    @property
    def bw(self):
        return self._bw

    @bw.setter
    def bw(self, x):
        self._bw[:, :] = x

    @property
    def bn(self):
        return self._bn

    @bn.setter
    def bn(self, x):
        self._bn[:, :] = x

    @property
    def bs(self):
        return self._bs

    @bs.setter
    def bs(self, x):
        self._bs[:, :] = x

#
# class scalarField(BaseField.baseField):
#
#     def __init__(self, data):
#         super().__init__(data)
#
#     @classmethod
#     def fromShape(cls, shape, value=0.0):
#         return cls( super().newField(shape,value) )
#
#     def __add__(self, other):
#         return scalarField( self.data + other.data )
#
#     def __mul__(self, other):
#         if isinstance(other, scalarField):
#             return scalarField(self.data * other.data)
#         elif isinstance(other, type(1.0)):
#             return scalarField(self.data * other)
#
#     def __sub__(self, other):
#         return scalarField(self.data - other.data)
#
#     def __neg__(self):
#         return scalarField(-self.data)

# merge with scalarField?
class varScalarField(ScalarField.scalarField):

    def __init__(self, mesh, geometry, internalValue=None):
        super().__init__( super().newField( mesh.shapeSimpleScalarField ) )
        self._boundary = dict.fromkeys(geometry.getBoundaryNames(), None)

    def setBoundaryCondition(self, boundaryName, boundaryType):
        self._boundary[boundaryName] = boundaryType

class vertexField:

    def __init__(self, mesh, value=0):
        self._mesh = mesh
        self._raw = PrimitiveFields.newVertexField(self._mesh, value)

        self._e = self._raw[:, 1:]
        self._w = self._raw[:, :-1]
        self._n = self._raw[:-1, :]
        self._s = self._raw[1:, :]

    @property
    def e(self):
        return self._e
    @property
    def w(self):
        return self._w
    @property
    def n(self):
        return self._n
    @property
    def s(self):
        return self._s




class varVectorField(vectorField):

    def __init__(self, mesh, geometry, value=0):

        super().__init__( *vectorField.UVfromMesh(mesh, value) )
        self._boundary = dict.fromkeys(geometry.getBoundaryNames(), None)

    def setBoundaryCondition(self, boundaryName, boundaryType, kwargs=None):
        self._boundary[boundaryName] = boundaryType
    #
    # def getSourceField(self):
    #     return self._b



def drawField(field, mesh):
    fieldType = type(field)
    if fieldType == varScalarField or fieldType == ScalarField.scalarField:
        drawCellField(field, mesh)
    elif fieldType == varVectorField or fieldType == vectorField:
        drawFaceField(field, mesh)

#
# def drawFaceField(field, mesh):
#
#     u_cell_primitive = Interpolation.getCellInterpolation(field._u, 'x')
#     u_cell = scalarField(u_cell_primitive)
#     #u_cell = scalarField(mesh=mesh, primitiveField=u_cell_primitive)
#     drawCellField(u_cell, mesh)
#
#     v_cell_primitive = Interpolation.getCellInterpolation(field._v, 'y')
#     v_cell = scalarField(v_cell_primitive)
# #    v_cell = scalarField(mesh=mesh, primitiveField=v_cell_primitive)
#     drawCellField(v_cell, mesh)

# a cfd solution method
def drawCellField(cellField, mesh):
    import matplotlib.pyplot as plt

    nbcellsX = mesh._cells_x
    nbcellsY = mesh._cells_y
    lenX = mesh._lenX
    lenY = mesh._lenY

    ax = plt.gca()
    map = ax.imshow(cellField._raw, cmap='hot', interpolation='nearest')

    ax.set_xticks(np.linspace(-0.5, nbcellsX-0.5, 5))
    ax.set_xticklabels( np.linspace(0,lenX,5))

    ax.set_yticks(np.linspace(-0.5, nbcellsY-0.5, 5))
    ax.set_yticklabels( np.linspace(0,lenY,5))

    plt.colorbar(map)
    plt.show()

