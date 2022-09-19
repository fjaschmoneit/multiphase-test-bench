import numpy as np
import importlib
#import scipy.sparse as sparse
#from scipy.sparse import dia_matrix, bsr_array
#from scipy.sparse.linalg import spsolve

import Fields
importlib.reload(Fields)

import ScalarField
importlib.reload(ScalarField)

class linearSystem:

    # constructor for cell linear system
    # consider to make constructor singleton
    def __init__(self, mesh, type):

        if type == 'scalar':
            self.nbCV = mesh._nbCells
            self.shape = (mesh._cells_y, mesh._cells_x)
        if type == 'vector_U':
            self.nbCV = mesh._nbCells-1
            self.shape = (mesh._cells_y, mesh._cells_x-1)

        # use scipy sprse amtrix classes here
        # self._A = dia_matrix( (nbCells, nbCells) ).toarray()
        # self._b = bsr_array( (nbCells)).toarray()

        self._A = np.zeros(self.shape)
        self._b = np.zeros((self.nbCV,))


    def solve(self):
        x = np.linalg.solve(self._A, self._b)
        return np.reshape(x, self.shape )


    # remove phi, that's only for fixing BCs, which I will move out of here
    # should not depend on mesh, should not depend on phi
    def update(self,mesh,F,D,S, phi):
    ### updates coefficient matrix and b vector from concatenated flux vectors
        self._A.fill(0.0)
        self._b.fill(0.0)

    # directional coefficient matrices are cellFields, very inefficient
    # here my difference schemes come into play
    # do I still need my fields here or can I just use primitive containers?
    # it is only due to the boundary conditions that I need scalarFields
        a_e = -ScalarField.scalarField(D.e - 0.5 * F.e)
        a_w = -ScalarField.scalarField(D.w + 0.5 * F.w)
        a_n = -ScalarField.scalarField(D.n - 0.5 * F.n)
        a_s = -ScalarField.scalarField(D.s + 0.5 * F.s)

        s_u = S
        a_p = -(a_e + a_w + a_n + a_s)

        # can I correct the boundary conditions when defining the input fluxes?
        # fixing boundary conditions
        if phi._boundary['top'] == 'zeroGradient':
            a_p.bn += a_n.bn
            a_n.bn = 0
        else:
            Tn = phi._boundary['top']
            s_u.bn -= a_n.bn * Tn
            a_n.bn = 0

        if phi._boundary['bottom'] == 'zeroGradient':
            a_p.bs += a_s.bs
            a_s.bs = 0
        else:
            Ts = phi._boundary['bottom']
            s_u.bs -= a_s.bs * Ts
            a_s.bs = 0

        if phi._boundary['left'] == 'zeroGradient':
            a_p.bw += a_w.bw
            a_w.bw = 0
        else:
            Tw = phi._boundary['left']
            s_u.bw -= a_w.bw * Tw
            a_w.bw = 0 #needed??

        if phi._boundary['right'] == 'zeroGradient':
            a_p.be += a_e.be
            a_e.be = 0
        else:
            Te = phi._boundary['right']
            s_u.be -= a_e.be * Te
            a_e.be = 0

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
        self._b = np.reshape(s_u._raw, s_u._raw.size)

        ny, nx = a_e.shape
        self._A = np.diag(a_p_serial) + np.diag(a_e_serial[:-1], 1) + np.diag(a_w_serial[1:], -1) + np.diag(a_s_serial[:-nx], nx) + np.diag(a_n_serial[nx:], -nx)
