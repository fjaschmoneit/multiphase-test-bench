import PrimitiveFields


def drawCellField(cellField):
    import matplotlib.pyplot as plt

    map = plt.imshow(cellField.raw, cmap='hot', interpolation='nearest')
    plt.colorbar(map)
    plt.show()


class parameterCellField():
    # a parameterCellField does not have a boundary
    def __init__(self, mesh, value=0, primitiveField=None):
        self.mesh = mesh
        self.type = 'parameterCellField'
        self.nx = mesh.cells_x
        self.ny = mesh.cells_y

        if primitiveField is not None:
            self.raw = primitiveField
        else:
            self.raw = PrimitiveFields.newCellField( mesh=mesh, value=value )

        self._bw = self.raw[:, :1]
        self._be = self.raw[:, -1:]
        self._bn = self.raw[:1, :]
        self._bs = self.raw[-1:, :]

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
        return parameterCellField( mesh=self.mesh, primitiveField=self.raw + other.raw )

    def __sub__(self, other):
        return parameterCellField(mesh=self.mesh, primitiveField=self.raw - other.raw)

    def __neg__(self):
        return parameterCellField(mesh=self.mesh, primitiveField=-self.raw)


    def fillWithConsecutiveValues(self):
        PrimitiveFields.fillWithConsecutiveValues(self.raw)

    def fillWithRandomIntegers(self):
        PrimitiveFields.fillWithRandomIntegers(self.raw)


class variableCellField(parameterCellField):

    def __init__(self, mesh, value=0, primitiveField=None):
        super().__init__(mesh, value, primitiveField)

        self._boundary = {}


class parameterFaceField:
    # access functions to underlying primitive structure
    def __init__(self, mesh, value=0):
        self.type = 'parameterFaceField'

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
