# mesh class should only hold meta data
import importlib

import Fields
#importlib.reload(Fields)

class cartesian2D():

    def __init__(self, len_x, len_y, res):
        ### a regular, rectangular mesh with nb of cells in [x,y] direction, according to dimensions parameter

        self._lenX = len_x
        self._lenY = len_y
        self._uniformSpacing = 1.0/res   # this is a scalar, which is only valid for regular cartesian meshes. reciprocal distances could be a member of every mesh
        self._cells_x = int(res*self._lenX)
        self._cells_y = int(res*self._lenY)
        self._nbCells = self._cells_x*self._cells_y

        self._invCellDist = None

        #--------------- public variables:
        self.shapeCVField = (self._cells_y, self._cells_x)
        self.shapeVarScalarField = (self._cells_y+1, self._cells_x+1)  # including a ghost frame
        self.shapeStaggered_U = (self._cells_y, self._cells_x+1)
        self.shapeStaggered_V = (self._cells_y+1, self._cells_x)
        #self.shapeFaceGrid = (self._cells_y+1, self._cells_x+1)

        #--------------- initializing variables:
        self.defineReciprocalDistances()


    def defineReciprocalDistances(self):
        ### sets recCellDist as a parameterFaceField with inverse cell distances between internal cells and face-cell distance at boundary

        ff = Fields.vectorField.fromShapes( self.shapeStaggered_U, self.shapeStaggered_V )
        ff.u.data =  1.0/self._uniformSpacing
        ff.v.data =  1.0/self._uniformSpacing

        ff.u.be *= 2
        ff.u.bw *= 2
        ff.v.bn *= 2
        ff.v.bs *= 2

        self._invCellDist = ff

    # def getCellVolumes(self):
    #     return Fields.scalarField(mesh=self, value=self._uniformSpacing**2)

    def getInverseCellDistances(self):
        return self._invCellDist

    def getFaceAreas(self):
        return Fields.vectorField.fromMesh(self, value=1.0)

    def getStats(self):
        print( self._nbCells )