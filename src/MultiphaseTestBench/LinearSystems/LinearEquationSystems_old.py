import numpy as np
#from scipy.sparse import diags

#from scipy.sparse import dia_matrix, bsr_array
#from scipy.sparse.linalg import spsolve

class linearSystem:

    def __init__(self, shape):
        # these dimensions correspond to a scalar control volume porblem.
        # staggered mesh problems are smaller, since only internal faces are solved for
        # in this case, the lin equation system is filled with zeros, but it doesn't change size
        self.shape = shape

        # just dummy matrix
        self.A = np.zeros((5,5))
        self.b = np.zeros((5,))

    def solve(self):
        # print(self.A)
        # print(self.b)
        x = np.linalg.solve(self.A, self.b)
        return np.reshape( x, self.shape )

    def reset(self):
        self.A.fill(0.0)
        self.b.fill(0.0)

    def set_b(self, b):
        linLength = self.shape[0] * self.shape[1]
        self.b = np.reshape(b, linLength)

    def set_e_coeffs(self,a):
        a_serialized = np.reshape(a, a.shape[0]*a.shape[1])
        self.A -= np.diag(a_serialized[:-1], 1)

    def set_w_coeffs(self, a):
        a_serialized = np.reshape(a, a.shape[0] * a.shape[1])
        self.A -= np.diag(a_serialized[1:], -1)

    def set_s_coeffs(self, a):
        ny,nx = a.shape
        a_serialized = np.reshape(a, nx*ny)
        self.A -= np.diag(a_serialized[:-nx], nx)

    def set_n_coeffs(self, a):
        ny, nx = a.shape
        a_serialized = np.reshape(a, nx * ny)
        self.A -= np.diag(a_serialized[nx:], -nx)

    def set_p_coeffs(self, a):
        a_serialized = np.reshape(a, a.shape[0] * a.shape[1])
        self.A += np.diag(a_serialized)

    # def update(self,F,D,Sc,Sp):
    # ### updates coefficient matrix and b vector from concatenated flux vectors
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
    #     linLength = self.shape[0]*self.shape[1]
    #
    #     self.b = np.reshape(Sc.data, linLength)
    #     return a_p
