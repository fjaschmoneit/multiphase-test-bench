# mesh class should only hold meta data
import Fields

class cartesian2D():

    def __init__(self, len_x=1, len_y=1, res = 10):
        ### a regular, rectangular mesh with nb of cells in [x,y] direction, according to dimensions parameter

        self.uniformSpacing = 1/res   # this is a scalar, which is only valid for regular cartesian meshes. reciprocal distances could be a member of every mesh
        self.dimensions = 2
        self.cells_x = int(res*len_x)
        self.cells_y = int(res*len_y)
        self.nbCells = self.cells_x*self.cells_y

        self.invCellDist = None
        self.defineReciprocalDistances()

    def defineReciprocalDistances(self):
        ### sets recCellDist as a parameterFaceField with inverse cell distances between internal cells and face-cell distance at boundary

        ff = Fields.parameterFaceField(mesh=self, value=1/self.uniformSpacing)

        ff.be = 2 / (self.uniformSpacing)
        ff.bw = 2 / (self.uniformSpacing)
        ff.bn = 2 / (self.uniformSpacing)
        ff.bs = 2 / (self.uniformSpacing)

        self.invCellDist = ff

