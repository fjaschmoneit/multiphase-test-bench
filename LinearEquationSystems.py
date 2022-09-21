import numpy as np

#import scipy.sparse as sparse
#from scipy.sparse import dia_matrix, bsr_array
#from scipy.sparse.linalg import spsolve

import Fields    # i keep this include until i have fixed te boundary conditions

class linearSystem:

    # constructor for cell linear system
    # consider to make constructor singleton
    def __init__(self, mesh):

        # use scipy sprse amtrix classes here
        # self._A = dia_matrix( (nbCells, nbCells) ).toarray()
        # self._b = bsr_array( (nbCells)).toarray()


        # these dimensions correspond to a scalar control volume porblem.
        # staggered mesh problems are smaller, since only internal faces are solved for
        # in this case, the lin equation system is filled with zeros, but it doesn't change size
        self.nbCV = mesh._nbCells
        self.shape = (mesh._cells_y, mesh._cells_x)

        self._A = np.zeros(self.shape)
        self._b = np.zeros((self.nbCV,))


    def solve(self):
        check why my matrix looks so funny

        x = np.linalg.solve(self._A, self._b)
        return np.reshape(x, self.shape )


    # remove phi, that's only for fixing BCs, which I will move out of here
    # should not depend on mesh, should not depend on phi
    def update(self,F,D,S, phi, fGov):
    ### updates coefficient matrix and b vector from concatenated flux vectors
        self._A.fill(0.0)
        self._b.fill(0.0)

    # directional coefficient matrices are cellFields, very inefficient
    # here my difference schemes come into play
    # do I still need my fields here or can I just use primitive containers?
    # it is only due to the boundary conditions that I need scalarFields

        a_e = -fGov.newScalarField(data=D.u.east - 0.5 * F.u.east)
        a_w = -fGov.newScalarField(data=D.u.west + 0.5 * F.u.west)
        a_n = -fGov.newScalarField(data=D.v.north - 0.5 * F.v.north)
        a_s = -fGov.newScalarField(data=D.v.south + 0.5 * F.v.south)

        s_c = S.Sc
        s_p = S.Sp
        a_p = -(a_e + a_w + a_n + a_s + s_p)

        # # can I correct the boundary conditions when defining the input fluxes?
        # # fixing boundary conditions
        # if phi._boundary['top'] == 'zeroGradient':
        #     a_p.bn += a_n.bn
        #     a_n.bn = 0
        # else:
        #     Tn = phi._boundary['top']
        #     s_c.bn -= a_n.bn * Tn
        #     a_n.bn = 0
        #
        # if phi._boundary['bottom'] == 'zeroGradient':
        #     a_p.bs += a_s.bs
        #     a_s.bs = 0
        # else:
        #     Ts = phi._boundary['bottom']
        #     s_c.bs -= a_s.bs * Ts
        #     a_s.bs = 0
        #
        # if phi._boundary['left'] == 'zeroGradient':
        #     a_p.bw += a_w.bw
        #     a_w.bw = 0
        # else:
        #     Tw = phi._boundary['left']
        #     s_c.bw -= a_w.bw * Tw
        #     a_w.bw = 0 #needed??
        #
        # if phi._boundary['right'] == 'zeroGradient':
        #     a_p.be += a_e.be
        #     a_e.be = 0
        # else:
        #     Te = phi._boundary['right']
        #     s_c.be -= a_e.be * Te
        #     a_e.be = 0

        # redefinition of my variables. should not be neccessary, when I only use primitive containers here
        a_e = a_e.data
        a_w = a_w.data
        a_n = a_n.data
        a_s = a_s.data
        a_p = a_p.data

        a_e_serial = np.reshape(a_e, a_e.size)
        a_w_serial = np.reshape(a_w, a_w.size)
        a_n_serial = np.reshape(a_n, a_n.size)
        a_s_serial = np.reshape(a_s, a_s.size)
        a_p_serial = np.reshape(a_p, a_p.size)

        # why is this reshaped?
        self._b = np.reshape(s_c._raw, s_c._raw.size)

        ny, nx = a_e.shape
        self._A = np.diag(a_p_serial) + np.diag(a_e_serial[:-1], 1) + np.diag(a_w_serial[1:], -1) + np.diag(a_s_serial[:-nx], nx) + np.diag(a_n_serial[nx:], -nx)
