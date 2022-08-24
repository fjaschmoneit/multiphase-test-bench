import Mesh, Fields, FlowModels

class Geometry:
    def __init__(self, lengthX, lengthY):
        self._lenX = lengthX
        self._lenY = lengthY
        self._boundaries = ['left', 'right', 'top', 'bottom']
        self._regions = ['internal']

    def getBoundaryNames(self):
        return self._boundaries

def createGeometry( typeName, kwargs ):
    if typeName == 'rectangle':
        return Geometry(*kwargs)

def createMesh(Geometry, res):
    len_x = Geometry._lenX
    len_y = Geometry._lenY
    return Mesh.cartesian2D(len_x = len_x, len_y= len_y, res=res)

class Simulation:
    def __init__(self, flowmodels, mesh, geometry):
        self._isCompiled = False
        self._mesh = mesh
        self._flowmodels = flowmodels.values()
        self._geometry = geometry

        self._variableFields = {}
        self._coefficientFields = {}
        self._scalarCoefficients = {}

        self.createfields()

    def createfields(self):
        self.deleteFields()
        for flowmodel in self._flowmodels:
            for var, type in flowmodel._variables.items():
                self._variableFields[var] = Fields.newField(type, self._mesh, self._geometry)

            for par, type in flowmodel._parameters.items():
                if par not in flowmodel._variables.keys():
                    self._coefficientFields[par] = Fields.newField(type, self._mesh)
        #
        # for v in self._variableFields:
        #     self._coefficientFields.pop(v, None)

    def initializeFieldData(self, mesh):
        for field in self._variableFields:
            field.initialize(mesh)

    def deleteFields(self):
        self._variableFields.clear()
        self._coefficientFields.clear()

    def getFields(self):
        return {**self._variableFields, **self._coefficientFields}

    def showfields(self):
        print("variable fields :")
        [print("\t", key, "\t", field._type, "\t", field._boundary) for key,field in self._variableFields.items()]
        print("parameter fields :")
        for key, field in self._coefficientFields.items():
            try:
                print("\t", key, "\t", field._type)
            except:
                print("\t", key, "\t", "scalar")
        print("\n")


    def updateLinSystems(self):
        # creating fields and their corresponding matrix equations

        # creates the linear equation systems
        for flowmodel in self._flowmodels:
            flowmodel.updateLinearEquationSystems(mesh=self._mesh, fields={**self._variableFields, **self._coefficientFields})

        self._isCompiled = True

    def display(self, field):
        Fields.drawField(field)



