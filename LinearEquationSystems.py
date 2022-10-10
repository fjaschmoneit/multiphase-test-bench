import numpy as np
#import scipy.sparse as sparse
#from scipy.sparse import dia_matrix, bsr_array
#from scipy.sparse.linalg import spsolve

class linearSystem:

    def __init__(self, mesh):

        # these dimensions correspond to a scalar control volume porblem.
        # staggered mesh problems are smaller, since only internal faces are solved for
        # in this case, the lin equation system is filled with zeros, but it doesn't change size
        self.nbCV = mesh._nbCells
        self.shape = None

        self._A = np.zeros((self.nbCV,self.nbCV))
        self._b = np.zeros((self.nbCV,))

    def solve(self):
        x = np.linalg.solve(self._A, self._b)
        return np.reshape( x, self.shape )

    def updatePressureEq(self, P, Sc):
        self._A.fill(0.0)
        # self._b.fill(0.0)      this is coupled with my source vector S, causing problems

        # implementation of differencing schemes here, currently only central difference
        a_e = -P.u.east
        a_w = -P.u.west
        a_n = -P.v.north
        a_s = -P.v.south

        a_p = -(a_e + a_w + a_n + a_s)

        self.shape = a_e.shape
        linLength = self.shape[0] * self.shape[1]  # I don't understand why I get this warning
        a_e_serial = np.reshape(a_e, linLength)
        a_w_serial = np.reshape(a_w, linLength)
        a_n_serial = np.reshape(a_n, linLength)
        a_s_serial = np.reshape(a_s, linLength)
        a_p_serial = np.reshape(a_p, linLength)

        self._b = np.reshape(Sc.data, linLength)

        ny, nx = a_e.shape
        self._A = np.diag(a_p_serial) + np.diag(a_e_serial[:-1], 1) + np.diag(a_w_serial[1:], -1) + np.diag(
            a_s_serial[:-nx], nx) + np.diag(a_n_serial[nx:], -nx)



    def set_e_coeffs(self,a_e):
        self.shape = a_e.shape
        linLength = self.shape[0]*self.shape[1]        # I don't understand why I get this warning
        a_e_serial = np.reshape(a_e, linLength)
        #ny, nx = a_e.shape
        #self._A += np.diag(a_e_serial[:-1], 1)
        self._A += np.diag(a_e_serial, 1)

    # not used?
    def reset(self):
        print('moin')
        self._A.fill(0.0)

    def update(self,F,D,Sc,Sp):
    ### updates coefficient matrix and b vector from concatenated flux vectors
        self._A.fill(0.0)
        #self._b.fill(0.0)      this is coupled with my source vector S, causing problems

        # implementation of differencing schemes here, currently only central difference
        a_e = -(D.u.east - 0.5 * F.u.east)
        a_w = -(D.u.west + 0.5 * F.u.west)
        a_n = -(D.v.north - 0.5 * F.v.north)
        a_s = -(D.v.south + 0.5 * F.v.south)

        a_p = -(a_e + a_w + a_n + a_s - Sp.data)

        self.shape = a_e.shape
        linLength = self.shape[0]*self.shape[1]        # I don't understand why I get this warning
        a_e_serial = np.reshape(a_e, linLength)
        a_w_serial = np.reshape(a_w, linLength)
        a_n_serial = np.reshape(a_n, linLength)
        a_s_serial = np.reshape(a_s, linLength)
        a_p_serial = np.reshape(a_p, linLength)

        self._b = np.reshape(Sc.data, linLength)

        ny, nx = a_e.shape
        self._A = np.diag(a_p_serial) + np.diag(a_e_serial[:-1], 1) + np.diag(a_w_serial[1:], -1) + np.diag(a_s_serial[:-nx], nx) + np.diag(a_n_serial[nx:], -nx)
