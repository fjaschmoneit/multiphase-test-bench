import LinearEquationSystems
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
        self._mesh = mesh
        self._geometry = geometry
        self._fieldRegistry = {}

        # this dict relates every field to its governing flowmodel
        self._fieldFlowModelLink = {}

        # should I not just have one? they are approx same size anyway
        self._scalarLinEqSystem = LinearEquationSystems.linearSystem(self._mesh, type='scalar')
        self._vectorLinEqSystem = LinearEquationSystems.linearSystem(self._mesh, type='vector_U')

        self._flowmodels = [fm(self._mesh, self._geometry, self._fieldRegistry, self._fieldFlowModelLink) for fm in flowmodels]

    def getFieldRegistry(self):
        return self._fieldRegistry

    def solveField(self, fieldname):

        flowmodel = self._fieldFlowModelLink[fieldname]
        flowmodel.updateLinearEquationSystem(self._mesh, self._fieldRegistry, self._scalarLinEqSystem)

        # flowModel dep
        self._fieldRegistry[fieldname]._raw[:,:] = self._scalarLinEqSystem.solve()
        #self._fieldRegistry[fieldname].internalEntriesEW = self._scalarLinEqSystem.solve()

    def display(self, field):
        Fields.drawField(field)



