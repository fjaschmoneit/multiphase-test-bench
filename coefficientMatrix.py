# nbCells = self.mesh.dimensions[0]
# self.field = np.zeros(nbCells)
self.coeffMatrix = np.zeros((mesh.cells_x, mesh.cells_y))  # not part of field
self.b = np.zeros(mesh.nbCells)
self.boundaries = dict()


def solve(self):
    self.field = np.linalg.solve(self.coeffMatrix, self.b)
