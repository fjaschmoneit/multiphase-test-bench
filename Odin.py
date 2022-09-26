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
        self._eqSystem = LinearEquationSystems.linearSystem(self._mesh)
        for fm in flowmodels:
            fm.initialize(self._mesh, self._fieldRegistry, self._fieldFlowModelLink)
            
    def solve(self, fieldname):

        field = self._fieldRegistry[fieldname]
        flowmodel = self._fieldFlowModelLink[fieldname]

        self.update( flowmodel, field )
        field.internal = self._eqSystem.solve()


    def update(self, flowmodel, field):
        flowmodel.updateFluxes(self._fieldRegistry)
        flowmodel.updateSourceField(self._fieldRegistry)
        flowmodel.correctBCs(self._fieldRegistry)

        self._eqSystem.update( F=flowmodel._convFluxes,
                               D=flowmodel._diffFluxes,
                               S=flowmodel._sourceField)


#--------- could also be external functions:

    def getFieldRegistry(self):
        return self._fieldRegistry


    def display(self, field, mesh):
        fGov = self._fieldRegistry['governor']
        fGov.drawField(field, mesh)

