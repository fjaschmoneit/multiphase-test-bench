import numpy as np


def newField(**kwargs):
    if 'data' in kwargs:
        data = kwargs.get('data')
        govModel = kwargs.get('governingModel')
        return baseField(data, govModel)
    elif 'shape' in kwargs:
        shape = kwargs.get('shape')
        value = kwargs.get('value', 0.0)
        # ghostSwitch = kwargs.get('includeGhostNodes', False)
        govModel = kwargs.get('governingModel')
        return baseField.fromShape(shape, value, govModel)
    else:
        print("newField must be called with type or data")

def drawField(field, mesh, title):
    drawCellField(field, mesh, title)

def drawCellField(cellField, mesh, title):
    import matplotlib.pyplot as plt

    nbcellsX = mesh.cells_x
    nbcellsY = mesh.cells_y
    lenX = mesh.lenX
    lenY = mesh.lenY

    ax = plt.gca()

    if isinstance(cellField, baseField):
        fieldData = cellField.data
    else:
        fieldData = cellField

    map = ax.imshow(fieldData, cmap='hot', interpolation='nearest')

    ax.set_title(title)
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

    def __init__(self, data, governingModel=None):

        self.govModel = governingModel          # remove also this.
        self.data = data

    #------------------ constructors
    @classmethod
    def fromShape(cls, shape, value=0.0, ghostSwitch=False, governingModel=None):
        field = cls(data = newDataField(shape, value), governingModel=governingModel )
        return field

    def fill(self, value):
        self.data[:,:] = value

    def getShape(self):
        return self.data.shape


class fieldContainer:
    def __init__(self, u, v):
        # these are base fields:
        self.u = u
        self.v = v

        self.govModel = u.govModel
        self.boundary = {}
