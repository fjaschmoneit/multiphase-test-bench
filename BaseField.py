import numpy as np

class baseField:

    def __init__(self, data):

        self._raw = data
        self.shape = self._raw.shape
        self._boundary = {}

        self._internalValues_U = self._raw[:,1:-1]
        self._internalValues_V = self._raw[1:-1,:]

        # directional values:
        self._east = self._raw[:, 1:]
        self._west = self._raw[:, :-1]
        self._north = self._raw[:-1, :]
        self._south = self._raw[1:, :]

        # ghost values:
        self._bw = self._raw[:, :1]
        self._be = self._raw[:, -1:]
        self._bn = self._raw[:1, :]
        self._bs = self._raw[-1:, :]

        # # boundary values:
        # self._bw = self._raw[:, 1:2]
        # self._be = self._raw[:, -2:-1]
        # self._bn = self._raw[1:2, :]
        # self._bs = self._raw[-2:-1, :]


    def __add__(self, other):
        return baseField(self.data + other.data)

    def __mul__(self, other):
        if isinstance(other, baseField):
            return baseField(self.data * other.data)
        elif isinstance(other, type(1.0)):
            return baseField(self.data * other)

    def __sub__(self, other):
        return baseField(self.data - other.data)

    def __neg__(self):
        return baseField(-self.data)

    #------------------ constructors
    @classmethod
    def fromShape(cls, shape, value=0.0):
        return cls( cls.newField(shape, value) )

    @classmethod
    def copy(cls, other):
        return cls( other.data )

    #----------------- static methods
    @staticmethod
    def newField(shape, value=0.0):
        field = np.ndarray(shape=shape, dtype=float, order='C')
        field.fill(value)
        return field

    # Should be done by flowmodel
    def setConstSource(self, value):
        self._source = value

    # Should be done by flowmodel
    def setBoundaryCondition(self, boundaryName, boundaryType):
        self._boundary[boundaryName] = boundaryType

    @property
    def data(self):
        return self._raw
    @data.setter
    def data(self,x):
        self._raw[:,:] = x


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

    @property
    def east(self):
        return self._east
    @east.setter
    def east(self, x):
        self._east[:,:] = x

    @property
    def west(self):
        return self._west
    @west.setter
    def west(self, x):
        self._west[:,:] = x

    @property
    def north(self):
        return self._north
    @north.setter
    def north(self, x):
        self._north[:,:] = x

    @property
    def south(self):
        return self._south
    @south.setter
    def south(self, x):
        self._south[:,:] = x
