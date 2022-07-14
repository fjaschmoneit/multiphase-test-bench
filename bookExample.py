# from: H K Versteeg and W Malalasekera, An Introduction to Computational Fluid Dynamics THE FINITE VOLUME METHOD, Second Edition, 2007
# pp 118

import numpy as np



class mesh(object):

    def __init__(self):
        self.cellCrds = {
            "ca": 0,
            "cb": 1,
            "cc": 2
        }

        self.faceCells = {
            "fa" : ("ca", "cb"),
            "fb" : ("cb", "cb"),
            "fc" : ("cb", ""),
            "fd" : ("ca", "")
        }



class coeffMatrix(object):

    def __init__(self, dim):
        self.A = np.zeros((dim,dim))
        self.b = np.zeros(dim)

    def setCoeffVector(self, t1, t2):
        self.b[0] = 200*t1
        self.b[1] = 0
        self.b[2] = 0
        self.b[3] = 0
        self.b[4] = 200*t2

    def setCoeffsMatrix(self):
        self.A[0, 0] = 300
        self.A[0, 1] = -100
        self.A[1, 0] = -100
        self.A[1, 1] = 200
        self.A[1, 2] = -100
        self.A[2, 1] = -100
        self.A[2, 2] = 200
        self.A[2, 3] = -100
        self.A[3, 2] = -100
        self.A[3, 3] = 200
        self.A[3, 4] = -100
        self.A[4, 3] = -100
        self.A[4, 4] = 300