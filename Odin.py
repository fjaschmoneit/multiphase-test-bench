import LinearEquationSystems
import Mesh
import Fields
import TransportModels
import PressureModels
from fieldAccess import *

class Geometry:
    def __init__(self, lengthX, lengthY):
        self.lenX = lengthX
        self.lenY = lengthY
        self.boundaries = [('left', 'west'), ('right', 'east'), ('top', 'north'), ('bottom', 'south')]
        self.regions = ['internal']

    def getBoundaryNames(self):
        return self.boundaries

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

def defineBoundaryCondition(field, boundaryName, boundaryType, fieldReg, **kwargs):
    value = kwargs.get('value', None)
    field.boundary[boundaryName] = (boundaryType, value)
    if boundaryType == 'fixedValue':
        if boundaryName == 'right':
            field.data[boundary_east] = value
        elif boundaryName == 'left':
            field.data[boundary_west] = value
        elif boundaryName == 'top':
            field.data[boundary_north] = value
        elif boundaryName == 'bottom':
            field.data[boundary_south] = value

#
# def setGoverningTransportModel(self, model):
#     self._govModel = model
#
#
# def getGoverningTransportModel(field, fieldReg):
#     fieldKey = getKeyFromValue(fieldReg, field)
#     return fieldReg[fieldKey]['governingModel']
#     #return self._govModel


# should directly change the b vector in the lin eq system
def updateSource(field, value, mesh):
    field.govModel.setSourceField(value * mesh.uniformSpacing)
#    fieldReg[fieldname]['governingModel'].setSourceField(value * mesh._uniformSpacing)
    # self._sourceField_c = self._fc.newField(type='scalarCV', value=value * self._mesh._uniformSpacing)


def solve(field):
    field.govModel.updateFluxes()
    field.govModel.updateSourceField()     i have to ingnore this for the shear stress test. fix it
    # field.govModel.correctBCs()
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
            self.fieldRegistry[fieldName].boundary = { name[0] :None for name in self.geometry.boundaries }


            # #= {
            #     'field':transportModelInstance.getDepField(),
            #     'governingModel': transportModelInstance,
            #     'boundaryConditions':{ name[0] :None for name in self._geometry._boundaries }
            # }

        for fieldName, fieldType in passiveFields.items():
            self.fieldRegistry[fieldName]= self.fieldCreator.newField(type=fieldType)

        fieldnames = self.fieldRegistry.keys()
        for fieldname in fieldnames:
            govModel = self.fieldRegistry[fieldname].govModel
            if govModel is not None:
                govModel.linkOtherFields(self.fieldRegistry)

            # if 'governingModel' in self._fieldRegistry[fieldname]:
            #     govModel = self._fieldRegistry[fieldname]['governingModel']
            #     govModel.linkOtherFields(self._fieldRegistry)
        # for transportmodel in flowmodels.values():
        #     transportmodel.linkOtherFields(self._fieldRegistry)



    #--------- could also be external functions:

    # def getFieldRegistry(self):
    #     return self.fieldRegistry

    def display(self, field, mesh):
        fGov = self.fieldCreator
        #fGov = self._fieldRegistrAy['governor']
        fGov.drawField(field, mesh)     # should not be a field creator method

