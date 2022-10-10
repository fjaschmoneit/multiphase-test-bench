import LinearEquationSystems
import Mesh
import Fields
import TransportModels
import PressureModels

class Geometry:
    def __init__(self, lengthX, lengthY):
        self._lenX = lengthX
        self._lenY = lengthY
        self._boundaries = [('left', 'west'), ('right', 'east'), ('top', 'north'), ('bottom', 'south')]
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
    def __init__(self, flowmodels, mesh, geometry, passiveFields={}):
        self._mesh = mesh
        self._geometry = geometry
        self._fieldRegistry = {}
        self._fieldCreator = Fields.fieldCreator(self._mesh)

        self._mesh.defineReciprocalDistances(self._fieldCreator, self._fieldRegistry)
        # this dict relates every field to its governing flowmodel
        self._fieldFlowModelLink = {}
        self._eqSystem = LinearEquationSystems.linearSystem(self._mesh)

        for fieldName, transportModelName in flowmodels.items():
            transportModelInstance = transportModelName(mesh, self._fieldCreator, self._fieldRegistry, self._eqSystem)
            flowmodels[fieldName] = transportModelInstance
            self._fieldRegistry[fieldName] = transportModelInstance.getDepField()

        for fieldName, fieldType in passiveFields.items():
            self._fieldRegistry[fieldName] = self._fieldCreator.newField(type=fieldType)

        for transportmodel in flowmodels.values():
            transportmodel.linkOtherFields(self._fieldRegistry)

            # fm(self._fieldCreator)
            # self._fieldFlowModelLink[]
            #fm.initializeFlowModel(self._mesh, self._fieldRegistry, self._fieldFlowModelLink)
            
    # def solve(self, fieldname):
    #
    #     field = self._fieldRegistry[fieldname]
    #     flowmodel = self._fieldFlowModelLink[fieldname]
    #
    #     self.update( flowmodel)
    #     field.data = self._eqSystem.solve()

        #sol = self._eqSystem.solve()
        # if isinstance(flowmodel, FlowModels.ScalarConvectionDiffusion):
        #     field.data = sol
        # elif isinstance(flowmodel, FlowModels.IncompressibleMomentumComp):
        #     field.data = sol

    # def update(self, flowmodel):
    #     #
    #     # flowmodel.updateFluxes(self._fieldRegistry)
    #     # flowmodel.updateSourceField(self._fieldRegistry, self._mesh)
    #     # flowmodel.correctBCs(self._fieldRegistry)
    #     #
    #     # Pressure system only needs one flux field:
    #     if isinstance(flowmodel, FlowModels.Pressure):
    #         flowmodel.updateFluxes_pressure(self._fieldRegistry, )
    #         flowmodel.updateSourceField(self._fieldRegistry, self._mesh)
    #         flowmodel.correctBCs(self._fieldRegistry)
    #         self._eqSystem.updatePressureEq( P=flowmodel._presFluxes, Sc=flowmodel._sourceField_c )
    #     else:
    #         flowmodel.updateFluxes(self._fieldRegistry)
    #         flowmodel.updateSourceField(self._fieldRegistry, self._mesh)
    #         flowmodel.correctBCs(self._fieldRegistry)
    #         self._eqSystem.update( F=flowmodel._convFluxes,
    #                            D=flowmodel._diffFluxes,
    #                            Sc=flowmodel._sourceField_c,
    #                            Sp=flowmodel._sourceField_p)

#--------- could also be external functions:

    def getFieldRegistry(self):
        return self._fieldRegistry

    def display(self, field, mesh):
        fGov = self._fieldCreator
        #fGov = self._fieldRegistrAy['governor']
        fGov.drawField(field, mesh)     # should not be a field creator method

