import numpy as np

class scalarField:

    def __init__(self, data):

        self._raw = data
        self._shape = self._raw.shape

        # internal values:
        self._e = self._raw[:, 1:]
        self._w = self._raw[:, :-1]
        self._n = self._raw[:-1, :]
        self._s = self._raw[1:, :]

        # boundary values:
        self._bw = self._raw[:, :1]
        self._be = self._raw[:, -1:]
        self._bn = self._raw[:1, :]
        self._bs = self._raw[-1:, :]
    #
    # class scalarField(BaseField.baseField):
    #
    #     def __init__(self, data):
    #         super().__init__(data)

    # @classmethod
    # def fromShape(cls, shape, value=0.0):
    #     return cls(newField(shape, value))

    def __add__(self, other):
        return scalarField(self.data + other.data)

    def __mul__(self, other):
        if isinstance(other, scalarField):
            return scalarField(self.data * other.data)
        elif isinstance(other, type(1.0)):
            return scalarField(self.data * other)

    def __sub__(self, other):
        return scalarField(self.data - other.data)

    def __neg__(self):
        return scalarField(-self.data)

    #------------------ constructors
    @classmethod
    def fromShape(cls, shape, value=0):
        return cls( cls.newField(shape, value) )

    @classmethod
    def copy(cls, other):
        return cls( other.data )

    #----------------- static methods
    @staticmethod
    def newField(shape, value=0):
        field = np.ndarray(shape=shape, dtype=float, order='C')
        field.fill(value)
        return field
    #
    # def fillWithConst(self, const):
    #     self._raw[:,:] = const
    #
    # def fillWithConsecutiveValues(self):
    #     PrimitiveFields.fillWithConsecutiveValues(self._raw)
    #
    # def fillWithRandomIntegers(self):
    #     PrimitiveFields.fillWithRandomIntegers(self._raw)

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
