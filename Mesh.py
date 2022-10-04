# mesh class should only hold meta data
import Fields

class cartesian2D():

    def __init__(self, len_x, len_y, res):
        ### a regular, rectangular mesh with nb of cells in [x,y] direction, according to dimensions parameter

        self._lenX = len_x
        self._lenY = len_y
        self._uniformSpacing = 1.0/res   # this is a scalar, which is only valid for regular cartesian meshes. reciprocal distances could be a member of every mesh
        self._cells_x = int(res*self._lenX)
        self._cells_y = int(res*self._lenY)
        self._nbCells = self._cells_x*self._cells_y


    def defineReciprocalDistances(self, fieldReg):
        ### sets recCellDist as a parameterFaceField with inverse cell distances between internal cells and face-cell distance at boundary

        fGov = fieldReg['governor']
        rCellDist = Fields.fieldContainer(
            u = fGov.newField(type='faces_u', value=1.0/self._uniformSpacing ),
            v = fGov.newField(type='faces_v', value=1.0/self._uniformSpacing )
        )
        fieldReg['invCellDist'] = rCellDist


    # def getCellVolumes(self):
    #     return Fields.scalarField(mesh=self, value=self._uniformSpacing**2)
    #
    # def getInverseCellDistances(self):
    #     return self._invCellDist

    def calcFaceAreas(self, fGov):
        fa = Fields.fieldContainer(
            u=fGov.newField(type='faces_u', value=1),
            v = fGov.newField(type='faces_v', value=1)
        )
        return fa

    def getStats(self):
        print( self._nbCells )