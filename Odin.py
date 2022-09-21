import LinearEquationSystems
import Mesh
from Fields import fieldGovernor
import FlowModels

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
        self._fieldRegistry = {
            'governor' : fieldGovernor(self._mesh)
        }

        self._mesh.defineReciprocalDistances(self._fieldRegistry)

        # this dict relates every field to its governing flowmodel
        self._fieldFlowModelLink = {}

        # why is the linear equation system defined here?
        # should I not just have one? they are approx same size anyway
        self._eqSystem = LinearEquationSystems.linearSystem(self._mesh)
        #self._vectorLinEqSystem = LinearEquationSystems.linearSystem(self._mesh, type='vector_U')

        #self._flowmodels = [fm(self._mesh, self._geometry, self._fieldRegistry, self._fieldFlowModelLink) for fm in flowmodels]
        for fm in flowmodels:
            fm.initializeFields(self._fieldRegistry, self._mesh)
            fm.initializeFluxesAndSource(self._fieldRegistry)
            fm.linkDepFieldToModel(self._fieldFlowModelLink)

    def solve(self, fieldname):

        field = self._fieldRegistry[fieldname]
        flowmodel = self._fieldFlowModelLink[fieldname]

        self.update( flowmodel, field )
        #flowmodel.updateLinearEquationSystem(self._mesh, self._fieldRegistry, self._scalarLinEqSystem)

        field.internal = self._eqSystem.solve()
        #self._fieldRegistry[fieldname].internalEntriesEW = self._scalarLinEqSystem.solve()


    def update(self, flowmodel, field):
        flowmodel.updateFluxes(self._fieldRegistry)
        flowmodel.updateSourceField(self._mesh, self._fieldRegistry)
#        S = flowmodel.calcSourceField(self._mesh, self._fieldRegistry)
        flowmodel.correctBCs()

        S = flowmodel._sourceField
        F = flowmodel._convFluxes
        D = flowmodel._diffFluxes
        self._eqSystem.update( F,D,S,field, self._fieldRegistry['governor'] )  # field is only passed because I don't treat the BCs in fluxes



#--------- could also be external functions:

    def getFieldRegistry(self):
        return self._fieldRegistry


    def display(self, field, mesh):
        fGov = self._fieldRegistry['governor']
        fGov.drawField(field, mesh)

