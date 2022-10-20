import numpy as np

class fieldCreator:

    def __init__(self, mesh):

        self._typeShapeDict = {
            'scalarCV'      : (mesh.cells_y, mesh.cells_x),
            'faces_u'       : (mesh.cells_y, mesh.cells_x + 1),
            'faces_v'       : (mesh.cells_y+1, mesh.cells_x),
            'staggered_u'   : (mesh.cells_y, mesh.cells_x+2),
            'vertex'        : (mesh.cells_y+1, mesh.cells_x+1),
        }

    def getFluxShapesFromCVShape(self, shape):
        [nb_v, nb_u] = list(shape)
        return ( (nb_v, nb_u+1), (nb_v+1, nb_u) )



    def newField(self, **kwargs):
        if 'data' in kwargs:
            data = kwargs.get('data')
            govModel = kwargs.get('governingModel')
            return baseField(data, govModel)
        elif 'shape' in kwargs:
            shape = kwargs.get('shape')
            value = kwargs.get('value', 0.0)
            ghostSwitch = kwargs.get('includeGhostNodes', False)
            govModel = kwargs.get('governingModel')
            return baseField.fromShape(shape, value, ghostSwitch, govModel)
        elif 'type' in kwargs:
            type = kwargs.get('type')
            shape = self._typeShapeDict[type]
            value = kwargs.get('value', 0.0)
            ghostSwitch = kwargs.get('includeGhostNodes', False)
            govModel = kwargs.get('governingModel')
            return baseField.fromShape(shape, value, ghostSwitch, govModel)
        else:
            print("newField must be called with type or data")

    @staticmethod
    def drawField(field, mesh):
        # fieldType = type(field)
        fieldCreator.drawCellField(field, mesh)
        # if fieldType == varScalarField or fieldType == ScalarField.scalarField:
        #     drawCellField(field, mesh)
        # elif fieldType == varVectorField or fieldType == vectorField:
        #     drawFaceField(field, mesh)

    @staticmethod
    def drawCellField(cellField, mesh):
        import matplotlib.pyplot as plt

        nbcellsX = mesh.cells_x
        nbcellsY = mesh.cells_y
        lenX = mesh.lenX
        lenY = mesh.lenY

        ax = plt.gca()
        map = ax.imshow(cellField.data, cmap='hot', interpolation='nearest')

        ax.set_xticks(np.linspace(-0.5, nbcellsX - 0.5, 5))
        ax.set_xticklabels(np.linspace(0, lenX, 5))

        ax.set_yticks(np.linspace(-0.5, nbcellsY - 0.5, 5))
        ax.set_yticklabels(np.linspace(0, lenY, 5))

        plt.colorbar(map)
        plt.show()


# ----------------- static methods
def newDataField(shape, value=0.0):
    field = np.ndarray(shape=shape, dtype=float, order='C')
    field.fill(value)
    return field

class baseField:

    def __init__(self, data, ghostSwitch=False, governingModel=None):

        self.govModel = governingModel
        self.data = data
        self.boundary = {}

    #------------------ constructors
    @classmethod
    def fromShape(cls, shape, value=0.0, ghostSwitch=False, governingModel=None):
        field = cls(data = newDataField(shape, value), ghostSwitch=ghostSwitch, governingModel=governingModel )
        return field

    def fill(self, value):
        self.data[:,:] = value



class fieldContainer:
    def __init__(self, u, v):
        # these are base fields:
        self.u = u
        self.v = v

        self.govModel = u.govModel
        self.boundary = {}
