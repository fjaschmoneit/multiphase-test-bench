import numpy as np
import scipy.sparse as sp
#from numpy import linalg as LA

#from scipy.sparse import diags
#from scipy.sparse import dia_matrix, bsr_array
#from scipy.sparse.linalg import spsolve

class linearSystem:

    def __init__(self, mesh):
        # these dimensions correspond to a scalar control volume porblem.
        # staggered mesh problems are smaller, since only internal faces are solved for
        # in this case, the lin equation system is filled with zeros, but it doesn't change size
        self.shape = None

        # just dummy matrix
        self.A = np.zeros((5,5))
        self.b = np.zeros((5,))

    def solve(self):

        self.filterTinyValues()

        # analytic matrix solving
        # x = np.linalg.solve(self.A, self.b)

        # sparse matrix solving
        A,b = self.convertToSparseMatrix()
        x = sp.linalg.spsolve(A, b)

        #print("solution ", x1.all() == x.all())
        return np.reshape( x, self.shape )

    def convertToSparseMatrix(self):
        # converting numpy nd array to sparse scipy matrix

        ny, nx = self.shape
        sprsA = sp.diags([self.A.diagonal(i) for i in [-nx, -1, 0, 1, nx]], [-nx, -1, 0, 1, nx]).tocsr()
        sprsb = sp.csr_matrix(self.b.reshape(-1, 1))

        # test if conversion is successful
        # print("coefficient matrix ", self.A.all() ==  sprsA.toarray().all())
        # print("b vector ", self.b.all() ==  sprsb.toarray().all())

        return sprsA,sprsb

    def filterTinyValues(self):
        lowThres = 1e-15
        self.A = np.where(np.absolute(self.A) > lowThres, self.A, 0)
        self.b = np.where(np.absolute(self.b) > lowThres, self.b, 0)


    def reset(self, shape):
        self.shape = shape
        nc= shape[0]*shape[1]
        self.A = np.zeros( (nc,nc) )
        self.b = np.zeros( nc )

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
    #
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
