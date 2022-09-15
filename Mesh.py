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

        self._invCellDist = self.defineReciprocalDistances()

    def defineReciprocalDistances(self):
        ### sets recCellDist as a parameterFaceField with inverse cell distances between internal cells and face-cell distance at boundary
        ff = Fields.vectorField(mesh=self, value=1/self._uniformSpacing)

        ff.be *= 2
        ff.bw *= 2
        ff.bn *= 2
        ff.bs *= 2

        return ff

    def getCellVolumes(self):
        return Fields.scalarField(mesh=self, value=self._uniformSpacing**2)

    def getInverseCellDistances(self):
        return self._invCellDist

    def getFaceAreas(self):
        #return Fields.vectorField(mesh=self, value=self._uniformSpacing**2)
        return Fields.vectorField(mesh=self, value=1)

    def getStats(self):
        print( self._nbCells )