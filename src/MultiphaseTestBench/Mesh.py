# mesh class should only hold meta data
import Fields
import MeshConfig


class cartesian2D():

    def __init__(self, len_x, len_y, res):
        ### a regular, rectangular mesh with nb of cells in [x,y] direction, according to dimensions parameter

        self.lenX = len_x
        self.lenY = len_y
        self.uniformSpacing = 1.0/res   # this is a scalar, which is only valid for regular cartesian meshes. reciprocal distances could be a member of every mesh
        self.cells_x = int(res*self.lenX)
        self.cells_y = int(res*self.lenY)
        self.nbCells = self.cells_x*self.cells_y

        MeshConfig.SHAPE_SCALAR_CV =            (self.cells_y, self.cells_x)
#        MeshConfig.SHAPE_SCALAR_CV_GHOST =      (self.cells_y+2, self.cells_x+2)
        MeshConfig.SHAPE_FACES_U =              (self.cells_y , self.cells_x + 1)
        MeshConfig.SHAPE_FACES_V =              (self.cells_y + 1, self.cells_x)
        MeshConfig.SHAPE_VERTEX =               (self.cells_y + 1, self.cells_x + 1)


    def calcInvCellDistance(self,direction):
        rDist = 1.0/self.uniformSpacing
        if direction == 'east' or direction == 'west':
            return Fields.newDataField(shape=MeshConfig.SHAPE_FACES_U, value=rDist)
        elif direction == 'north' or direction == 'south':
            return Fields.newDataField(shape=MeshConfig.SHAPE_FACES_V, value=rDist)


    def calcInvCellDistances(self):
        ### sets recCellDist as a parameterFaceField with inverse cell distances between internal cells and face-cell distance at boundary

        f_u = self.calcInvCellDistance('west')
        f_v = self.calcInvCellDistance('south')
        return (f_u, f_v)


    def calcFaceArea(self, direction):
        constArea = self.uniformSpacing**2
#        constArea = 1
        if direction == 'east' or direction == 'west':
            return Fields.newDataField(shape=MeshConfig.SHAPE_FACES_U, value=constArea)
        elif direction == 'north' or direction == 'south':
            return Fields.newDataField(shape=MeshConfig.SHAPE_FACES_V, value=constArea)

    def calcFaceAreas(self):
        # f_u = Fields.newDataField( shape=fGov.typeShapeDict['faces_u'], value=1.0 )
        # f_v = Fields.newDataField( shape=fGov.typeShapeDict['faces_v'], value=1.0 )
        f_u = self.calcFaceArea('east')
        f_v = self.calcFaceArea('south')

        return (f_u, f_v)

    def getStats(self):
        print( self.nbCells )