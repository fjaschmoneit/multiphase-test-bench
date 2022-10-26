import LinearEquationSystems
import Mesh
import Fields
import TransportModels
import PressureModels
from fieldAccess import *
import numpy as np

class Geometry:

    boundaries = {
        'left': 'west',
        'right': 'east',
        'top': 'north',
        'bottom': 'south'}

    def __init__(self, lengthX, lengthY):
        self.lenX = lengthX
        self.lenY = lengthY
        self.regions = ['internal']

    def getBoundaryNames(self):
        return self.boundaries

    @ classmethod
    def getCompDirectionFromName(cls, name):
        return cls.boundaries[name]

def createGeometry( typeName, kwargs ):
    if typeName == 'rectangle':
        return Geometry(*kwargs)

def createMesh(Geometry, res):
    len_x = Geometry.lenX
    len_y = Geometry.lenY
    return Mesh.cartesian2D(len_x = len_x, len_y= len_y, res=res)

def getKeyFromValue(dict,value):
    for k,v in dict.items():
        if v == value:
            return k
    print("error:\t{value} not fount in dictionary {dict}".format(**locals()))

def defineBoundaryCondition(field, boundaryName, **argDict):

    # linking BC in transport model:
    dir = Geometry.getCompDirectionFromName(boundaryName)
    argDict['direction'] = dir
    field.boundary[boundaryName] =  argDict

    # setting initial values (only possible, when data is initialized)
    (boundary_dir, boundary_nb1_dir) = fieldSlice( dir )
    boundaryType = argDict.get('type')
    if boundaryType == 'fixedValue':
        # in scalar transport, setting face BC values on cell centres
        field.data[boundary_dir] = argDict.get('value')
    elif boundaryType == 'zeroGradient':
        neighborField = field.data[boundary_nb1_dir]
        if neighborField.size != 0:         # this could happen in 1D simulations
            field.data[boundary_dir] = neighborField

def calcCollocatedVelocityField(u,v):
    col_u = 0.5*(u[west] + u[east])
    col_v = 0.5*(v[north] + v[south])

    #return np.dstack((col_u,col_v))
    return (col_u, col_v)

def calcVelocityMagnitude( vel ):
    #col_u, col_v = np.dsplit(vel, 2)
    col_u = vel[0]
    col_v = vel[1]
    return np.sqrt( col_u**2 + col_v**2 )


def listAvailableBoundaryModels(field):
    return field.govModel.listAvailableBoundaryModels()

# should directly change the b vector in the lin eq system
def setConstSource(field, value, mesh):
    field.govModel.setConstSourceField(  Fields.newDataField(field.getShape(), value*mesh.uniformSpacing ) )

def solve(field):
    field.govModel.updateFluxes()
    field.govModel.updateSourceField()
    field.govModel.updateLinSystem()
    return field.govModel.solve()

class Simulation:
    def __init__(self, flowmodels, mesh, geometry, passiveFields={}):
        self.mesh = mesh
        self.geometry = geometry
        self.fieldRegistry = {}
        self.fieldCreator = Fields.fieldCreator(self.mesh)

        self.mesh.defineReciprocalDistances(self.fieldCreator, self.fieldRegistry)
        self._eqSystem = LinearEquationSystems.linearSystem(self.mesh)

        for fieldName, transportModelName in flowmodels.items():
            transportModelInstance = transportModelName(mesh, self.fieldCreator, self.fieldRegistry, self._eqSystem)

            self.fieldRegistry[fieldName] = transportModelInstance.getDepField()
            self.fieldRegistry[fieldName].govModel = transportModelInstance
#            self.fieldRegistry[fieldName].boundary = { name[0] :None for name in self.geometry.boundaries }
            self.fieldRegistry[fieldName].boundary = { name :None for name in self.geometry.boundaries }

        for fieldName, fieldType in passiveFields.items():
            self.fieldRegistry[fieldName]= self.fieldCreator.newField(type=fieldType)

        fieldnames = self.fieldRegistry.keys()
        for fieldname in fieldnames:

            if isinstance( self.fieldRegistry[fieldname], Fields.baseField):
                govModel = self.fieldRegistry[fieldname].govModel
                if govModel is not None:
                    govModel.linkOtherFields(self.fieldRegistry)

    def display(self, field, mesh):
        fGov = self.fieldCreator
        fGov.drawField(field, mesh)

