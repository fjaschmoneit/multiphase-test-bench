import numpy as np
import BaseField


class fieldGovernor:

    def __init__(self, mesh):
        self.shapeCVField = (mesh._cells_y, mesh._cells_x)
        self.shapeVarScalarField = (mesh._cells_y + 2, mesh._cells_x + 2)  # including a ghost frame
        self.shapeFaces_u = (mesh._cells_y, mesh._cells_x + 1)
        self.shapeFaces_v = (mesh._cells_y + 1, mesh._cells_x)

    @classmethod
    def sourceField(cls, fieldGovernor):
        return sourceFields(fieldGovernor)

    @classmethod
    def newScalarField(cls, **kwargs):
        if 'data' in kwargs:
            return scalarField(kwargs['data'])
        elif 'shape' in kwargs:
            if 'value' not in kwargs:
                value = 0.0
            else:
                value = kwargs['value']
            return scalarField.fromShape(kwargs['shape'], value)
        else:
            print("error: incorrect constructor arguments supplied.")

    @classmethod
    def newVectorField(cls,**kwargs):
        if 'u' in kwargs and 'v' in kwargs:
            return vectorField(kwargs['u'], kwargs['v'])
        elif 'shape_u' in kwargs and 'shape_v' in kwargs:
            if 'value' not in kwargs:
                value = 0.0
            else:
                value = kwargs['value']
            u = BaseField.baseField.fromShape(shape=kwargs['shape_u'], value=value)
            v = BaseField.baseField.fromShape(shape=kwargs['shape_v'], value=value)

            # u = cls.newScalarField(shape=kwargs['shape_u'], value=value)
            # v = cls.newScalarField(shape=kwargs['shape_v'], value=value)
            return vectorField(u,v)
        else:
            print("error: incorrect constructor arguments supplied.")


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
        map = ax.imshow(cellField._raw, cmap='hot', interpolation='nearest')

        ax.set_xticks(np.linspace(-0.5, nbcellsX - 0.5, 5))
        ax.set_xticklabels(np.linspace(0, lenX, 5))

        ax.set_yticks(np.linspace(-0.5, nbcellsY - 0.5, 5))
        ax.set_yticklabels(np.linspace(0, lenY, 5))

        plt.colorbar(map)
        plt.show()

class sourceFields:
    def __init__(self, fGov):
        # make these fGov methods
        self.Sc = BaseField.baseField.fromShape(shape=fGov.shapeCVField, value=0.0)
        self.Sp = BaseField.baseField.fromShape(shape=fGov.shapeCVField, value=0.0)

    def setConstantSource(self, value):
        self.Sc.data = value

class scalarField(BaseField.baseField):
    # a scalar field has ghost nodes
    def __init__(self, data):

        super().__init__(data)

        self.internal = self._raw[1:-1,1:-1]

        # ghost values:
        self._gw = self._raw[:, :1]
        self._ge = self._raw[:, -1:]
        self._gn = self._raw[:1, :]
        self._gs = self._raw[-1:, :]

        # boundary values:
        self._bw = self._raw[:, 1:2]
        self._be = self._raw[:, -2:-1]
        self._bn = self._raw[1:2, :]
        self._bs = self._raw[-2:-1, :]

    @classmethod
    def fromShape(cls, shape, value=0.0):
        data = BaseField.baseField.newField(shape, value)
        return cls( data )

    @property
    def ge(self):
        return self._ge
    @ge.setter
    def ge(self,x):
        self._ge[:,:] = x

    @property
    def gw(self):
        return self._gw
    @gw.setter
    def gw(self, x):
        self._gw[:, :] = x

    @property
    def gn(self):
        return self._gn
    @gn.setter
    def gn(self, x):
        self._gn[:, :] = x

    @property
    def gs(self):
        return self._gs
    @gs.setter
    def gs(self, x):
        self._gs[:, :] = x

class vectorField:
    def __init__(self, u, v):
        # these are scalar fields:
        self._u = u
        self._v = v

    @classmethod
    def fromShapes(cls, shape_u, shape_v, value=0.0):
        u = BaseField.baseField.fromShape(shape_u, value)
        v = BaseField.baseField.fromShape(shape_v, value)
        return cls(u, v)

#-------- defining algebra
    def __add__(self, other):
        u = self._u + other._u
        v = self._v + other._v
        return vectorField(u,v)

    def __mul__(self, other):
        if isinstance(other, vectorField):
            return vectorField(self._u * other._u, self._v * other._v)
        elif isinstance(other, type(1.0)):
            return vectorField(self._u * other, self._v * other)
        else:
            print("multiplication not defined")

    #-------- defining all getters and setters:

    @property
    def u(self):
        return self._u
    @u.setter
    def u(self, x):
        self._u.data[:, :] = x
    @property
    def v(self):
        return self._v
    @v.setter
    def v(self, x):
        self._v.data[:, :] = x
    #
    # @property
    # def be(self):
    #     return self._be
    # @be.setter
    # def be(self, x):
    #     self._be = x
    #
    # @property
    # def bw(self):
    #     return self._bw
    # @bw.setter
    # def bw(self, x):
    #     self._bw = x
    #
    # @property
    # def bn(self):
    #     return self._bn
    # @bn.setter
    # def bn(self, x):
    #     self._bn = x
    #
    # @property
    # def bs(self):
    #     return self._bs
    # @bs.setter
    # def bs(self, x):
    #     self._bs = x


# merge with vector field
class varVectorField(vectorField):

    def __init__(self, mesh, geometry, value=0):

        super().__init__( *vectorField.UVfromMesh(mesh, value) )
        self._boundary = dict.fromkeys(geometry.getBoundaryNames(), None)

    def setBoundaryCondition(self, boundaryName, boundaryType, kwargs=None):
        self._boundary[boundaryName] = boundaryType
    #
    # def getSourceField(self):
    #     return self._b



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
