


class baseField:

    def __init__(self, rawField):

        self._raw = rawField
        self._shape = self._raw.shape

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

