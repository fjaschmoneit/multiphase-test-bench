# mesh class should only hold meta data
import Fields

class cartesian2D():

    def __init__(self, len_x, len_y, res):
        ### a regular, rectangular mesh with nb of cells in [x,y] direction, according to dimensions parameter

        self.lenX = len_x
        self.lenY = len_y
        self.uniformSpacing = 1.0/res   # this is a scalar, which is only valid for regular cartesian meshes. reciprocal distances could be a member of every mesh
        self.cells_x = int(res*self.lenX)
        self.cells_y = int(res*self.lenY)
        self.nbCells = self.cells_x*self.cells_y


    def defineReciprocalDistances(self, fGov, fieldReg):
        ### sets recCellDist as a parameterFaceField with inverse cell distances between internal cells and face-cell distance at boundary

        # #fGov = fieldReg['governor']
        # rCellDist = Fields.fieldContainer(
        #     u = fieldCreator.newField(type='faces_u', value=1.0/self.uniformSpacing ),
        #     v = fieldCreator.newField(type='faces_v', value=1.0/self.uniformSpacing )
        # )

        f_u = Fields.newDataField( shape=fGov.typeShapeDict['faces_u'], value=1.0/self.uniformSpacing )
        f_v = Fields.newDataField( shape=fGov.typeShapeDict['faces_v'], value=1.0/self.uniformSpacing )
        fieldReg['invCellDist'] = (f_u, f_v)


    # def getCellVolumes(self):
    #     return Fields.scalarField(mesh=self, value=self._uniformSpacing**2)
    #
    # def getInverseCellDistances(self):
    #     return self._invCellDist

    def calcFaceAreas(self, fGov):
        # u=fGov.newField(type='faces_u', value=1),
        # v = fGov.newField(type='faces_v', value=1)
        f_u = Fields.newDataField( shape=fGov.typeShapeDict['faces_u'], value=1.0 )
        f_v = Fields.newDataField( shape=fGov.typeShapeDict['faces_v'], value=1.0 )
        return (f_u, f_v)

    def getStats(self):
        print( self.nbCells )