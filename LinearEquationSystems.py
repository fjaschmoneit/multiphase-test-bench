import numpy as np
import Fields

# remove phi, that's only for fixing BCs, which I will move out of here
# should not depend on mesh
def createCoefficientMatrix(mesh,F,D,S, phi):
### creates coefficient matrix and b vector from concatenated flux vectors

    # directional coefficient matrices are cellFields, very inefficient
    # here my difference schemes come into play
    # do I still need my fields here or can I just use primitive containers?
    a_e = -Fields.parameterCellField(mesh=mesh, primitiveField=D.e - 0.5 * F.e)
    a_w = -Fields.parameterCellField(mesh=mesh, primitiveField=D.w + 0.5 * F.w)
    a_n = -Fields.parameterCellField(mesh=mesh, primitiveField=D.n - 0.5 * F.n)
    a_s = -Fields.parameterCellField(mesh=mesh, primitiveField=D.s + 0.5 * F.s)

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
        a_w.bw = 0

    if phi._boundary['right'] == 'zeroGradient':
        a_p.be += a_e.be
        a_e.be = 0
    else:
        Te = phi._boundary['right']
        s_u.be -= a_e.be * Te
        a_e.be = 0

    # redefinition of my variables. should not be neccessary, when I only use primitive containers here
    a_e = a_e._raw
    a_w = a_w._raw
    a_n = a_n._raw
    a_s = a_s._raw
    a_p = a_p._raw
    S = s_u._raw

    ny,nx = a_e.shape

    a_e_serial = np.reshape(a_e, a_e.size)
    a_w_serial = np.reshape(a_w, a_w.size)
    a_n_serial = np.reshape(a_n, a_n.size)
    a_s_serial = np.reshape(a_s, a_s.size)
    a_p_serial = np.reshape(a_p, a_p.size)
    S_serial = np.reshape(S, S.size)

    A = np.diag(a_p_serial)
    A = A+np.diag(a_e_serial[:-1], 1)
    A = A+np.diag(a_w_serial[1:], -1)
    A = A+np.diag(a_s_serial[:-nx], nx)
    A = A+np.diag(a_n_serial[nx:], -nx)

    return A,S_serial

def solveLinearSystem(A,b):
    return np.linalg.solve(A, b)
