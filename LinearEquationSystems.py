import numpy as np
from scipy.sparse import diags

#from scipy.sparse import dia_matrix, bsr_array
#from scipy.sparse.linalg import spsolve

class linearSystem:

    def __init__(self, mesh):
        # these dimensions correspond to a scalar control volume porblem.
        # staggered mesh problems are smaller, since only internal faces are solved for
        # in this case, the lin equation system is filled with zeros, but it doesn't change size
        self._shape = None

        # just dummy matrix
        self._A = np.zeros((5,5))
        self._b = np.zeros((5,))

    def solve(self):
        x = np.linalg.solve(self._A, self._b)
        return np.reshape( x, self._shape )

    def reset(self, shape):
        # self._A.fill(0.0)
        # self._b.fill(0.0)
        self._shape = shape
        nc= shape[0]*shape[1]
        self._A = np.zeros( (nc,nc) )
        self._b = np.zeros( nc )

    def set_b(self, b):
        linLength = self._shape[0] * self._shape[1]
        self._b = np.reshape(b, linLength)

    def set_e_coeffs(self,a):
        a_serialized = np.reshape(a, a.shape[0]*a.shape[1])
        self._A -= np.diag(a_serialized[:-1], 1)

    def set_w_coeffs(self, a):
        a_serialized = np.reshape(a, a.shape[0] * a.shape[1])
        self._A -= np.diag(a_serialized[1:], -1)

    def set_s_coeffs(self, a):
        ny,nx = a.shape
        a_serialized = np.reshape(a, nx*ny)
        self._A -= np.diag(a_serialized[:-nx], nx)

    def set_n_coeffs(self, a):
        ny, nx = a.shape
        a_serialized = np.reshape(a, nx * ny)
        self._A -= np.diag(a_serialized[nx:], -nx)

    def set_p_coeffs(self, a):
        a_serialized = np.reshape(a, a.shape[0] * a.shape[1])
        self._A += np.diag(a_serialized)

    def update(self,F,D,Sc,Sp):
    ### updates coefficient matrix and b vector from concatenated flux vectors

        # implementation of differencing schemes here, currently only central difference
        a_e = -(D.u.east - 0.5 * F.u.east)
        a_w = -(D.u.west + 0.5 * F.u.west)
        a_n = -(D.v.north - 0.5 * F.v.north)
        a_s = -(D.v.south + 0.5 * F.v.south)

        a_p = -(a_e + a_w + a_n + a_s - Sp.data)

        self.set_e_coeffs(a_e)
        self.set_w_coeffs(a_w)
        self.set_s_coeffs(a_s)
        self.set_n_coeffs(a_n)
        self.set_p_coeffs(a_p)

        linLength = self._shape[0]*self._shape[1]

        # make this also a funtion
        self._b = np.reshape(Sc.data, linLength)
        return a_p

    #
    # def updatePressure(self, P, Sc, Sp):
    #     ### updates coefficient matrix and b vector from concatenated flux vectors
    #
    #     # implementation of differencing schemes here, currently only central difference
    #     a_e = -(D.u.east - 0.5 * F.u.east)
    #     a_w = -(D.u.west + 0.5 * F.u.west)
    #     a_n = -(D.v.north - 0.5 * F.v.north)
    #     a_s = -(D.v.south + 0.5 * F.v.south)
    #
    #     a_p = -(a_e + a_w + a_n + a_s - Sp.data)
    #
    #     self.set_e_coeffs(a_e)
    #     self.set_w_coeffs(a_w)
    #     self.set_s_coeffs(a_s)
    #     self.set_n_coeffs(a_n)
    #     self.set_p_coeffs(a_p)
    #
    #     linLength = self._shape[0] * self._shape[1]
    #
    #     # make this also a funtion
    #     self._b = np.reshape(Sc.data, linLength)

    # def updatePressureEq(self, P, Sc):
    #     self._A.fill(0.0)
    #     # self._b.fill(0.0)      this is coupled with my source vector S, causing problems
    #
    #     # implementation of differencing schemes here, currently only central difference
    #     a_e = -P.u.east
    #     a_w = -P.u.west
    #     a_n = -P.v.north
    #     a_s = -P.v.south
    #
    #     a_p = -(a_e + a_w + a_n + a_s)
    #
    #     self.shape = a_e.shape
    #     linLength = self.shape[0] * self.shape[1]
    #     a_e_serial = np.reshape(a_e, linLength)
    #     a_w_serial = np.reshape(a_w, linLength)
    #     a_n_serial = np.reshape(a_n, linLength)
    #     a_s_serial = np.reshape(a_s, linLength)
    #     a_p_serial = np.reshape(a_p, linLength)
    #
    #     self._b = np.reshape(Sc.data, linLength)
    #
    #     ny, nx = a_e.shape
    #     self._A = np.diag(a_p_serial) + np.diag(a_e_serial[:-1], 1) + np.diag(a_w_serial[1:], -1) + np.diag(
    #         a_s_serial[:-nx], nx) + np.diag(a_n_serial[nx:], -nx)
