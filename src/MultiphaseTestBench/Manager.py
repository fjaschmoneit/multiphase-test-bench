import LinearEquationSystems
import Mesh
import Fields
import MeshConfig
from fieldAccess import *
import numpy as np
import ObjectRegistry as objReg
import TransportModels
import TransportModelsSWE
import PressureModels

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
    """
    creating a geometry object

    :param typeName: ['rectangle']
    :param kwargs: see Geometry constructor arguments
    :return: Geometry object
    """

    if typeName == 'rectangle':
        return Geometry(*kwargs)

def createMesh(Geometry, res):
    """
    creating a mesh from a geometry with specified resolution

    :param Geometry: geometry object
    :param res: integer [cells per meter]
    :return: Mesh object
    """

    len_x = Geometry.lenX
    len_y = Geometry.lenY
    return Mesh.cartesian2D(len_x = len_x, len_y= len_y, res=res)

def getKeyFromValue(dict,value):
    for k,v in dict.items():
        if v == value:
            return k
    print("error:\t{value} not fount in dictionary {dict}".format(**locals()))

def defineBoundaryCondition(field, boundaryName, type, **argDict):
    """
    defining a boundary conditions on a specified boundary and field

    :param type: ['zeroGradient', 'fixedValue']
    :param field: variable flied, e.g. U,p, ...
    :param boundaryName: name to boundary
    :param argDict: value for bcType=='fixedValue'
    :return: None
    """
    # linking BC in transport model:
    dir = Geometry.getCompDirectionFromName(boundaryName)
    argDict['type'] = type
    argDict['direction'] = dir
    field.govModel.boundary[boundaryName] =  argDict

    # # setting initial values (only possible, when data is initialized)
    # (boundary_dir, boundary_nb1_dir) = fieldSlice( dir )
    # boundaryType = argDict.get('type')
    # if boundaryType == 'fixedValue':
    #     # in scalar transport, setting face BC values on cell centres
    #     field.data[boundary_dir] = argDict.get('value')
    # elif boundaryType == 'zeroGradient':
    #     neighborField = field.data[boundary_nb1_dir]
    #     if neighborField.size != 0:         # this could happen in 1D simulations
    #         field.data[boundary_dir] = neighborField

def calcCollocatedVelocityField(u,v):
    """
    Interpolates a collocated velocity field from a staggered velocity field.

    :param u: staggered velocity component
    :param v: staggered velocity component
    :return: collocated velocity field
    """
    col_u = 0.5*(u[west] + u[east])
    col_v = 0.5*(v[north] + v[south])

    #return np.dstack((col_u,col_v))
    return (col_u, col_v)

def calcVelocityMagnitude( vel ):
    """
    Calculates the velocity magnitude.

    :param vel: collocated velocity field
    :return: cellField with velocity magnitudes
    """
    #col_u, col_v = np.dsplit(vel, 2)
    col_u = vel[0]
    col_v = vel[1]
    return np.sqrt( col_u**2 + col_v**2 )

def calcFlowRate(field, boundaryName):
    """
    Calculates the flow rate through a given boundary.

    :param field: field name
    :param boundaryName: name of boundary
    :return: net flow rate through boundary
    """
    direction = field.boundary[boundaryName]['direction']
    (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)

    vel = field.data[boundary_dir]
    fa = field.govModel.mesh.calcFaceArea(direction)[boundary_dir]

    return np.sum( vel*fa )


def listAvailableBoundaryModels(field):
    """
    Returns a list of available boundary condition models for the given field.

    :param field: field name
    :return: list of available boundary models
    """
    return field.govModel.listAvailableBoundaryModels()

# should directly change the b vector in the lin eq system
def setConstSource(field, value, mesh):
    """
    defining a constant scalar source in the entire domain.

    :param field: scalar field
    :param value: source per volume
    :param mesh: the computational mesh object
    :return: None
    """
    # shape = field.getShape()[internal].shape
    field.govModel.setConstSourceField(  value*mesh.uniformSpacing )

def transientSolve(field, dt, method):
    """
    solve the field as to the field's transient transport model

    :param field: scalar or vector field
    :return: None
    """
    field.govModel.updateFluxes()
    field.govModel.updateSourceField()
    field.govModel.updateTransientCoeffs(dt, method)


    field.govModel.updateLinSystem()
    return field.govModel.solve()



def solve(field):
    """
    solve the field as to the field's transport model

    :param field: scalar or vector field
    :return: None
    """
#    field.govModel.updateBoundaries()
    field.govModel.updateFluxes()
    field.govModel.updateSourceField()
    field.govModel.updateLinSystem()
    return field.govModel.solve()

def getField(name):
    """
    retrieve a field reference from object registry

    :param name: field name
    :return: reference to field object
    """

    return objReg.FIELDS[name]

def initialize(flowmodels, mesh, geometry, closure={}, passiveFields={}):
    """
    Initializing all flow models and includes corresponding fields to global object registry.

    :param flowmodels: dictionary of fields and corresponding flow models
    :param mesh: mesh object
    :param geometry: geometry object
    :param passiveFields: dictionary of additional fields and their respective field type
    :return: None
    """
    if closure is None:
        closure = {}
    objReg.MESH = mesh
    objReg.GEOM = geometry
    objReg.LINEAR_SYSTEM = LinearEquationSystems.linearSystem(objReg.MESH)

    # initializing flow models and corresponding fields:
    for fieldName, transportModelName in flowmodels.items():
        transportModelInstance = transportModelName(objReg.MESH, objReg.LINEAR_SYSTEM)

        objReg.FIELDS[fieldName] = transportModelInstance.getField()
        objReg.FIELDS[fieldName].govModel = transportModelInstance
        objReg.FIELDS[fieldName].govModel.boundary = {name: None for name in objReg.GEOM.boundaries}

    # initializing passive fields:
    for fieldName, fieldType in passiveFields.items():
        typeShapeDict = {
            'faces_u': MeshConfig.SHAPE_FACES_U,
            'faces_v': MeshConfig.SHAPE_FACES_V,
            'scalarCV': MeshConfig.SHAPE_SCALAR_CV
        }
        objReg.FIELDS[fieldName] = Fields.newField(shape=typeShapeDict[fieldType])

    # linking fields between flow models:
    fieldnames = objReg.FIELDS.keys()
    for fieldname in fieldnames:
        if isinstance(objReg.FIELDS[fieldname], Fields.baseField):
            govModel = objReg.FIELDS[fieldname].govModel
            if govModel is not None:
                if len(closure.keys())==0:
                    govModel.linkOtherFields([])
                else:
                    govModel.linkOtherFields(closure[fieldname])


def display(field, mesh, title):
    """
    Rudimental visual data analysis tool. Fields are plotted using matplotlib.

    :param title: title of figure
    :param field: field name
    :param mesh: mesh object
    :return: None
    """
    Fields.drawField(field, mesh, title)
