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
        self.shape = (mesh._cells_y, mesh._cells_x)

        self._A = np.zeros(self.shape)
        self._b = np.zeros((self.nbCV,))

    def solve(self):
        return np.reshape( np.linalg.solve(self._A, self._b), self.shape )

    def update(self,F,D,S):
    ### updates coefficient matrix and b vector from concatenated flux vectors
        self._A.fill(0.0)
        #self._b.fill(0.0)      this is coupled with my source vector S, causing problems

        # implementation of differencing schemes here, currently only central difference
        a_e = -(D.u.east - 0.5 * F.u.east)
        a_w = -(D.u.west + 0.5 * F.u.west)
        a_n = -(D.v.north - 0.5 * F.v.north)
        a_s = -(D.v.south + 0.5 * F.v.south)

        a_p = -(a_e + a_w + a_n + a_s - S.Sp.data)

        linLength = a_e.size        # I don't understand why I get this warning
        a_e_serial = np.reshape(a_e, linLength)
        a_w_serial = np.reshape(a_w, linLength)
        a_n_serial = np.reshape(a_n, linLength)
        a_s_serial = np.reshape(a_s, linLength)
        a_p_serial = np.reshape(a_p, linLength)

        self._b = np.reshape(S.Sc.data, linLength)

        ny, nx = a_e.shape
        self._A = np.diag(a_p_serial) + np.diag(a_e_serial[:-1], 1) + np.diag(a_w_serial[1:], -1) + np.diag(a_s_serial[:-nx], nx) + np.diag(a_n_serial[nx:], -nx)
