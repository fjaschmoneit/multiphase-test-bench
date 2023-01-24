import numpy as np

class linearSystem:

    def __init__(self, dim):
        self.shape = dim
        self.nby, self.nbx = self.shape
        self.nbtot = self.nbx*self.nby

        self.a = np.zeros((self.nbtot,self.nbtot))
        self.b = np.zeros((self.nby,self.nbx))

        self.diagonal = {
            'west' : (np.arange(self.nbtot - 1) + 1, np.arange(self.nbtot - 1)),
            'east' : (np.arange(self.nbtot - 1), np.arange(self.nbtot - 1) + 1),
            'north': (np.arange(self.nbtot - self.nbx) + self.nbx, np.arange(self.nbtot - self.nbx)),
            'south': (np.arange(self.nbtot - self.nbx), np.arange(self.nbtot - self.nbx) + self.nbx),
            'centre': (np.arange(self.nbtot), np.arange(self.nbtot))
        }

    def reset(self):
        self.a.fill(0.0)
        self.b.fill(0.0)

    def getCoefficientsInDirection(self,direction):
        return self.a[self.diagonal[direction]]

    def getCoefficientsInOppositeDirection(self, direction):
        if direction == 'west':
            opdirection = 'east'
        elif direction == 'east':
            opdirection = 'west'
        elif direction == 'north':
            opdirection = 'south'
        elif direction == 'south':
            opdirection = 'north'
        else:
            print("ERROR: unknown direction.")
            opdirection = 0

        return self.a[self.diagonal[opdirection]]


# coefficient getters:
    @property
    def a_w(self):
#        return self.a[self.diagonal['west']]
        return np.reshape(self.a[self.diagonal['west']], (self.nby, self.nbx-1))

    @property
    def a_e(self):
        return np.reshape(self.a[self.diagonal['east']], (self.nby, self.nbx-1))
#        return self.a[self.diagonal['east']]

    @property
    def a_n(self):
        return self.a[self.diagonal['north']]

    @property
    def a_s(self):
        return self.a[self.diagonal['south']]

    @property
    def a_p(self):
        return np.reshape(self.a[self.diagonal['centre']], (self.nby,self.nbx))
        #return self.a[self.diagonal['centre']]

# coefficient setters:
    @a_w.setter
    def a_w(self, dirFlux):
        if np.isscalar(dirFlux):
            self.a[self.diagonal['west']] = dirFlux
        else:
            self.a[self.diagonal['west']] = np.reshape(dirFlux, self.nbtot)[1:]
            #self.a[self.diagonal_w] = np.reshape(dirFlux, (self.nbx-1)*self.nby)
#            self.a[self.diagonal_w] = np.reshape(dirFlux, self.nbtot)[1:]

    @a_e.setter
    def a_e(self, dirFlux):
        if np.isscalar(dirFlux):
            self.a[self.diagonal['east']] = dirFlux
        else:
            self.a[self.diagonal['east']] = np.reshape(dirFlux, self.nbtot)[:-1]
#            self.a[self.diagonal_e] = np.reshape(dirFlux, self.nbtot)

    @a_n.setter
    def a_n(self, dirFlux):
        if np.isscalar(dirFlux):
            self.a[self.diagonal['north']] = dirFlux
        else:
           self.a[self.diagonal['north']] = np.reshape(dirFlux, self.nbtot)[self.nbx:]
#            self.a[self.diagonal_n] = np.reshape(dirFlux, self.nbtot)

    @a_s.setter
    def a_s(self, dirFlux):
        if np.isscalar(dirFlux):
            self.a[self.diagonal['south']] = dirFlux
        else:
            self.a[self.diagonal['south']] = np.reshape(dirFlux, self.nbtot)[:-self.nbx]
 #          self.a[self.diagonal_s] = np.reshape(dirFlux, self.nbtot)

    @a_p.setter
    def a_p(self, dirFlux):
        if np.isscalar(dirFlux):
            self.a[self.diagonal['centre']] = dirFlux
        else:
            self.a[self.diagonal['centre']] = np.reshape(dirFlux, self.nbtot)

    def solve(self):
        b = np.reshape(self.b, self.nbtot)
        x = np.linalg.solve(self.a, b)
        return np.reshape(x, self.shape)
