import numpy as np



# implement these free data views?
def east(data):
    return data[:, 1:]

def west(data):
    return data[:, :-1]

class baseField:

    def __init__(self, data, ghostSwitch=False):

        self._data = data
        self._internal = self._data[1:-1, 1:-1]

        if ghostSwitch:
            shapelist = list(self._data.shape)
            shapeWithGhost = tuple( [shapelist[0]+2, shapelist[1]+2] )
            self._raw = self.newDataField( shape=shapeWithGhost )
            self._raw[1:-1, 1:-1] = self._data

            # ghost access:
            self._gw = self._raw[:, :1]
            self._ge = self._raw[:, -1:]
            self._gn = self._raw[:1, :]
            self._gs = self._raw[-1:, :]
        else:
            self._raw = self._data

        self._boundary = {}

        self._internal_u = self._data[:,1:-1]
        self._internal_v = self._data[1:-1,:]

        # directional values:
        self._east = self._data[:, 1:]
        self._west = self._data[:, :-1]
        self._north = self._data[:-1, :]
        self._south = self._data[1:, :]

        # boundary values:
        self._bw = self._data[:, :1]
        self._be = self._data[:, -1:]
        self._bn = self._data[:1, :]
        self._bs = self._data[-1:, :]

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
    def fromShape(cls, shape, value=0.0, ghostSwitch=False):
        field = cls(data = cls.newDataField(shape, value), ghostSwitch=ghostSwitch )
        return field

    # ----------------- static methods
    @staticmethod
    def newDataField(shape, value=0.0):
        field = np.ndarray(shape=shape, dtype=float, order='C')
        field.fill(value)
        return field

    # a copy constructor only for the base class makes little sense
    @classmethod
    def copy(cls, other):
        return cls( other.data )

    def fill(self, value):
        self._data[:,:] = value

    def defineBoundaryCondition(self, boundaryName, boundaryType, **kwargs):
        value = kwargs.get('value', None)
        self._boundary[boundaryName] = (boundaryType, value)

        if boundaryType == 'fixedValue':
            if boundaryName == 'top':
                self.bn.fill( value )
            elif boundaryName == 'bottom':
                self.bs.fill( value )
            elif boundaryName == 'right':
                self.be.fill( value )
            elif boundaryName == 'left':
                self.bw.fill(value)


    # # Should be done by flowmodel
    # def setBoundaryCondition(self, boundaryName, boundaryType, value):
    # def setBoundaryCondition(self, kwargs):
    #     self._boundary[boundaryName] = boundaryType
    #     if boundaryType != 'zeroGradient':  # fix boundary conditions
    #         self.bw = boundaryType

    @property
    def internal(self):
        return self._internal
    @internal.setter
    def internal(self, x):
        self._internal[:, :] = x

    @property
    def raw(self):
        return self._raw
    @raw.setter
    def raw(self, x):
        self._raw[:, :] = x

    @property
    def data(self):
        return self._data
    @data.setter
    def data(self,x):
        self._data[:,:] = x

    @property
    def internal_u(self):
        return self._internal_u
    @internal_u.setter
    def internal_u(self,x):
        self._internal_u[:,:] = x

    @property
    def internal_v(self):
        return self._internal_v
    @internal_v.setter
    def internal_v(self, x):
        self._internal_v[:, :] = x

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
