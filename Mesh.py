# mesh class should only hold meta data
import Fields

class cartesian2D():

    def __init__(self, len_x, len_y, res):
        ### a regular, rectangular mesh with nb of cells in [x,y] direction, according to dimensions parameter

        self._lenX = len_x
        self._lenY = len_y
        self._uniformSpacing = 1/res   # this is a scalar, which is only valid for regular cartesian meshes. reciprocal distances could be a member of every mesh
        #self._dimensions = 2
        self._cells_x = int(res*self._lenX)
        self._cells_y = int(res*self._lenY)
        self._nbCells = self._cells_x*self._cells_y

        self.invCellDist = None
        self.defineReciprocalDistances()

    def defineReciprocalDistances(self):
        ### sets recCellDist as a parameterFaceField with inverse cell distances between internal cells and face-cell distance at boundary

        ff = Fields.parameterFaceField(mesh=self, value=1/self._uniformSpacing)

        ff.be = 2 / (self._uniformSpacing)
        ff.bw = 2 / (self._uniformSpacing)
        ff.bn = 2 / (self._uniformSpacing)
        ff.bs = 2 / (self._uniformSpacing)

        self.invCellDist = ff

    def getCellVolumes(self):
        return Fields.parameterCellField(mesh=self, value=self._uniformSpacing**3)

    def getStats(self):
        print( self._nbCells )