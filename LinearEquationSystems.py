import numpy as np


def createCoefficientMatrix(a_e,a_w,a_n,a_s,a_p, S):
### creates coefficient matrix and b vector from concatenated flux vectors

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
