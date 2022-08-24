import PrimitiveFields, LinearEquationSystems, Interpolation
import numpy as np


def drawField(field):
    fieldType = type(field)
    if fieldType == variableCellField or fieldType == parameterCellField:
        drawCellField(field)
    elif fieldType == variableFaceField or fieldType == parameterFaceField:
        drawFaceField(field)


def drawFaceField(field):
    mesh = field._mesh

    u_cell_primitive = Interpolation.getCellInterpolation(field._u, 'x')
    u_cell = parameterCellField(mesh=mesh, primitiveField=u_cell_primitive)
    drawCellField(u_cell)

    v_cell_primitive = Interpolation.getCellInterpolation(field._v, 'y')
    v_cell = parameterCellField(mesh=mesh, primitiveField=v_cell_primitive)
    drawCellField(v_cell)

# a cfd solution method
def drawCellField(cellField):
    import matplotlib.pyplot as plt

    nbcellsX = cellField._mesh._cells_x
    nbcellsY = cellField._mesh._cells_y
    lenX = cellField._mesh._lenX
    lenY = cellField._mesh._lenY

    ax = plt.gca()
    map = ax.imshow(cellField._raw, cmap='hot', interpolation='nearest')

    ax.set_xticks(np.linspace(-0.5, nbcellsX-0.5, 5))
    ax.set_xticklabels( np.linspace(0,lenX,5))

    ax.set_yticks(np.linspace(-0.5, nbcellsY-0.5, 5))
    ax.set_yticklabels( np.linspace(0,lenY,5))

    plt.colorbar(map)
    plt.show()


def newField(fieldType, mesh, geometry=None):
    if fieldType == 'vectorField':
        if geometry is None:
            return parameterFaceField(mesh)
#            return parameterFaceField()
        else:
            return variableFaceField(mesh, geometry)
#            return variableFaceField()
    elif fieldType == 'scalarField':
        if geometry is None:
            return parameterCellField(mesh)
            #return parameterCellField()
        else:
            return variableCellField(mesh, geometry)
#            return variableCellField()

    elif fieldType == 'scalar':
        return 0.0
    else:
        print("unknown field type: ", fieldType)


class parameterCellField():

    def __init__(self, mesh, value=0, primitiveField=None):
        self._mesh = mesh
        self._nx = mesh._cells_x
        self._ny = mesh._cells_y

        if primitiveField is not None:
            self._raw = primitiveField
        else:
            self._raw = PrimitiveFields.newCellField( mesh=mesh, value=value )

        self._bw = self._raw[:, :1]
        self._be = self._raw[:, -1:]
        self._bn = self._raw[:1, :]
        self._bs = self._raw[-1:, :]

    @property
    def be(self):
        return self._be
    @be.setter
    def be(self,x):
        self._be[:,:] = x

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

    def __add__(self, other):
        return parameterCellField( mesh=self._mesh, primitiveField=self._raw + other._raw )

    def __mul__(self, other):
        if isinstance(other, parameterCellField):
            return parameterCellField(mesh=self._mesh, primitiveField=self._raw * other._raw)
        elif isinstance(other, type(1.0)):
            return parameterCellField(mesh=self._mesh, primitiveField=self._raw * other)
    # def __truediv__(self, other):
    #     return parameterCellField(mesh=self._mesh, primitiveField=self._raw / other._raw)

    def __sub__(self, other):
        return parameterCellField(mesh=self._mesh, primitiveField=self._raw - other._raw)

    def __neg__(self):
        return parameterCellField(mesh=self._mesh, primitiveField=-self._raw)

    def fillWithConsecutiveValues(self):
        PrimitiveFields.fillWithConsecutiveValues(self._raw)

    def fillWithRandomIntegers(self):
        PrimitiveFields.fillWithRandomIntegers(self._raw)


class variableCellField(parameterCellField):

    def __init__(self, mesh, geometry, internalValue=None):
        super().__init__(mesh)
        self._type = 'variableCellField'
        self._boundary = dict.fromkeys(geometry.getBoundaryNames(), None)
        self._A = None
        self._b = None
        self.internalValue = internalValue

    def setBoundaryCondition(self, boundaryName, boundaryType, kwargs=None):
        self._boundary[boundaryName] = boundaryType

    def solve(self):
        x = LinearEquationSystems.solveLinearSystem(self._A, self._b)
        self._raw[:, :] = x.reshape(self._ny, self._nx)


class parameterFaceField:
    # access functions to underlying primitive structure
    def __init__(self, mesh, value=0, primitiveFieldX=None, primitiveFieldY=None):
        self._type = 'parameterFaceField'
        self._mesh = mesh

        if(primitiveFieldX is None and primitiveFieldY is None):
            self._u = PrimitiveFields.newFaceField_x( self._mesh, value )
            self._v = PrimitiveFields.newFaceField_y(self._mesh, value)
        else:
            self._u = primitiveFieldX
            self._v = primitiveFieldY

        self._u_internal = self._u[:,1:-1]
        self._v_internal = self._v[1:-1, :]

        self._e = self._u[:, 1:]
        self._w = self._u[:, :-1]
        self._n = self._v[:-1,:]
        self._s = self._v[1:,:]

        self._be = self._u[:, -1:]
        self._bw = self._u[:, :1]
        self._bn = self._v[:1, :]
        self._bs = self._v[-1:, :]

    def setInternalValue(self, valueVector):
        self._u[:,:] = valueVector[0]
        self._v[:,:] = valueVector[1]

#-------- defining algebra
    def __add__(self, other):
        primitiveFieldX = self._u + other._u
        primitiveFieldY = self._v + other._v
        return parameterFaceField(mesh=self._mesh, primitiveFieldX=primitiveFieldX, primitiveFieldY=primitiveFieldY )

    def __mul__(self, other):
        if isinstance(other, parameterFaceField):
            return parameterFaceField(mesh=self._mesh, primitiveFieldX=self._u * other._u, primitiveFieldY=self._v * other._v)
        elif isinstance(other, type(1.0)):
            return parameterFaceField(mesh=self._mesh, primitiveFieldX=self._u * other, primitiveFieldY=self._v * other)
        else:
            print("multiplication not defined")

    # def __truediv__(self, other):
    # #     return parameterCellField(mesh=self._mesh, primitiveField=self._raw / other._raw)
    #
    # def __sub__(self, other):
    #     return parameterCellField(mesh=self._mesh, primitiveField=self._raw - other._raw)
    #
    # def __neg__(self):
    #     return parameterCellField(mesh=self._mesh, primitiveField=-self._raw)

    #-------- defining all getters and setters:

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
        self._be[:,:] = x

    @property
    def bw(self):
        return self._bw
    @bw.setter
    def bw(self, x):
        self._bw[:,:] = x

    @property
    def bn(self):
        return self._bn
    @bn.setter
    def bn(self, x):
        self._bn[:,:] = x

    @property
    def bs(self):
        return self._bs
    @bs.setter
    def bs(self, x):
        self._bs[:,:] = x


class variableFaceField(parameterFaceField):

    def __init__(self, mesh, geometry, value=0, primitiveField=None):
        super().__init__(mesh, value, primitiveField)

        self._type = 'variableFaceField'
        self._boundary = dict.fromkeys(geometry.getBoundaryNames(), None)

    def setBoundaryCondition(self, boundaryName, boundaryType, kwargs=None):
        self._boundary[boundaryName] = boundaryType


if __name__ == '__main__':
    import Mesh

    mesh = Mesh.cartesian2D(0.5, 0.1, 10)
    a = parameterFaceField(mesh=mesh, value=2)
    b = parameterFaceField(mesh=mesh, value=2)
    c = a+b
    d = parameterFaceField(mesh=mesh, value=4)
    # i have to implement an equality operator
    print( c==d )