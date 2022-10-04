import numpy as np
import BaseField

class fieldGovernor:

    def __init__(self, mesh):
        # self.shapeCVField = (mesh._cells_y, mesh._cells_x)
        # self.shapeVarScalarField = (mesh._cells_y + 2, mesh._cells_x + 2)  # including a ghost frame
        # self.shapeFaces_u = (mesh._cells_y, mesh._cells_x + 1)
        # self.shapeFaces_v = (mesh._cells_y + 1, mesh._cells_x)
        # self.shapeVerticesInternal_u = (mesh._cells_y+1, mesh._cells_x-1)
        # self.shapeFacesInternal_u = (mesh._cells_y, mesh._cells_x-1)
        #
        # self.shapeVertices = (mesh._cells_y+1, mesh._cells_x+1)
        # # self.shapeStaggerdFaces_u = (mesh._cells_y, mesh._cells_x-1)
        # # self.shapeStaggerdFaces_v = (mesh._cells_y-1, mesh._cells_x)

        self._typeShapeDict = {
            'scalarCV'      : (mesh._cells_y, mesh._cells_x),
            'faces_u'       : (mesh._cells_y, mesh._cells_x + 1),
            'faces_v'       : (mesh._cells_y+1, mesh._cells_x),
            'internalVertices_u' : (mesh._cells_y+1, mesh._cells_x-1),
            'faceSource_u'  : (mesh._cells_y, mesh._cells_x-1),
            'faceSource_v'  : (mesh._cells_y-1, mesh._cells_x)
        }

    def newField(self, **kwargs):
        if 'data' in kwargs:
            data = kwargs.get('data')
            return BaseField.baseField(data)
        elif 'type' in kwargs:
            type = kwargs.get('type')
            shape = self._typeShapeDict[type]
            value = kwargs.get('value', 0.0)
            ghostSwitch = kwargs.get('includeGhostNodes', False)
            return BaseField.baseField.fromShape(shape, value, ghostSwitch)
        else:
            print("newField must be called with type or data")

    @staticmethod
    def drawField(field, mesh):
        # fieldType = type(field)
        fieldGovernor.drawCellField(field, mesh)
        # if fieldType == varScalarField or fieldType == ScalarField.scalarField:
        #     drawCellField(field, mesh)
        # elif fieldType == varVectorField or fieldType == vectorField:
        #     drawFaceField(field, mesh)

    @staticmethod
    def drawCellField(cellField, mesh):
        import matplotlib.pyplot as plt

        nbcellsX = mesh._cells_x
        nbcellsY = mesh._cells_y
        lenX = mesh._lenX
        lenY = mesh._lenY

        ax = plt.gca()
        map = ax.imshow(cellField.data, cmap='hot', interpolation='nearest')

        ax.set_xticks(np.linspace(-0.5, nbcellsX - 0.5, 5))
        ax.set_xticklabels(np.linspace(0, lenX, 5))

        ax.set_yticks(np.linspace(-0.5, nbcellsY - 0.5, 5))
        ax.set_yticklabels(np.linspace(0, lenY, 5))

        plt.colorbar(map)
        plt.show()




class fieldContainer:
    def __init__(self, u, v):
        # these are base fields:
        self._u = u
        self._v = v

#-------- defining algebra
    def __add__(self, other):
        u = self._u + other._u
        v = self._v + other._v
        return fieldContainer(u,v)

    def __mul__(self, other):
        if isinstance(other, fieldContainer):
            return fieldContainer(self._u * other._u, self._v * other._v)
        elif isinstance(other, type(1.0)):
            return fieldContainer(self._u * other, self._v * other)
        else:
            print("multiplication not defined")

    #-------- defining all getters and setters:

    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, x):
        self._u = x

    @property
    def v(self):
        return self._v
    @v.setter
    def v(self, x):
        self._v = x

#
# def drawFaceField(field, mesh):
#
#     u_cell_primitive = Interpolation.getCellInterpolation(field._u, 'x')
#     u_cell = scalarField(u_cell_primitive)
#     #u_cell = scalarField(mesh=mesh, primitiveField=u_cell_primitive)
#     drawCellField(u_cell, mesh)
#
#     v_cell_primitive = Interpolation.getCellInterpolation(field._v, 'y')
#     v_cell = scalarField(v_cell_primitive)
# #    v_cell = scalarField(mesh=mesh, primitiveField=v_cell_primitive)
#     drawCellField(v_cell, mesh)
