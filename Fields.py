import PrimitiveFields, LinearEquationSystems
import numpy as np

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
        else:
            return variableFaceField(mesh, geometry)
    elif fieldType == 'scalarField':
        if geometry is None:
            return parameterCellField(mesh)
        else:
            return variableCellField(mesh, geometry)

    elif fieldType == 'scalar':
        return 0.0
    else:
        print("unknown field type: ", fieldType)


class parameterCellField():
    # a parameterCellField does not have a boundary
    def __init__(self, mesh, value=0, primitiveField=None):
        self._mesh = mesh
        self._type = 'parameterCellField'
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

    def __init__(self, mesh, geometry, value=0, primitiveField=None):
        super().__init__(mesh, value, primitiveField)
        self._type = 'variableCellField'
        self._boundary = dict.fromkeys(geometry.getBoundaryNames(), None)
        self._A = None
        self._b = None


    def setBoundaryCondition(self, boundaryName, boundaryType, kwargs=None):
        self._boundary[boundaryName] = boundaryType

    def solve(self):
        x = LinearEquationSystems.solveLinearSystem(self._A, self._b)
        self._raw[:, :] = x.reshape(self._ny, self._nx)


class parameterFaceField:
    # access functions to underlying primitive structure
    def __init__(self, mesh, value=0, primitiveField=None):
        self._type = 'parameterFaceField'

        self._entries_EW = PrimitiveFields.newFaceField_x( mesh, value )
        self._entries_NS = PrimitiveFields.newFaceField_y(mesh, value)

        self._internalEntries_EW = self._entries_EW[:,1:-1]
        self._internalEntries_NS = self._entries_NS[1:-1, :]

        self._e = self._entries_EW[:, 1:]
        self._w = self._entries_EW[:, :-1]
        self._n = self._entries_NS[:-1,:]
        self._s = self._entries_NS[1:,:]

        self._be = self._entries_EW[:, -1:]
        self._bw = self._entries_EW[:, :1]
        self._bn = self._entries_NS[:1, :]
        self._bs = self._entries_NS[-1:, :]

#-------- defining all getters and setters:

    @property
    def entries_NS(self):
        return self._entries_NS

    @entries_NS.setter
    def entries_NS(self, x):
        self._entries_NS[:, :] = x

    @property
    def entries_EW(self):
        return self._entries_EW

    @entries_EW.setter
    def entries_EW(self, x):
        self._entries_EW[:, :] = x

    @property
    def internalEntries_NS(self):
        return self._internalEntries_NS
    @internalEntries_NS.setter
    def internalEntries_NS(self, x):
        self._internalEntries_NS[:, :] = x

    @property
    def internalEntries_EW(self):
        return self._internalEntries_EW
    @internalEntries_EW.setter
    def internalEntries_EW(self, x):
        self._internalEntries_EW[:, :] = x

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
