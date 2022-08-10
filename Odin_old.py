class Simulation:
    def __init__(self, geometry, mesh, flowmodels, kernel):
        self._geometry = geometry
        self._flowmodels = flowmodels
        self._kernel = kernel
        self._continuumproperties = Modules.ContinuumProperties(self._flowmodels, self._geometry)
        self._mesh = mesh

    def initialize(self):
        self._kernel.initialize(self._flowmodels, self._continuumproperties, self._mesh)

# should all these be kernel classes??
class Modules:
    # class Geometry:
    #     def __init__(self, lengthX, lengthY):
    #         self._lenX = lengthX
    #         self._lenY = lengthY
    #         self._boundaries = ['left', 'right', 'top', 'bottom']
    #         self._regions = ['internal']

    # class FlowModel:
    #     # ask kernel for available flow models
    #     class ScalarConvectionDiffusion():
    #         def __init__(self, scalarField='c', velocityField='U', diffusionCoefficient='D'):
    #             self._variables = {
    #                 scalarField: 'scalarField'
    #             }
    #             self._parameters = {
    #                 velocityField: 'vectorField',
    #                 diffusionCoefficient:'scalar'
    #             }
    #
    #         def show(self):
    #             print("ScalarConvectionDiffusion:\n\tvariables   {o._variables}\n\tparameters  {o._parameters}\n".format(o=self))
    #
    #     class IncompressibleFlow():
    #         def __init__(self, velocityField='U', pressureField='p', dynViscosity='nu'):
    #
    #             self._variables = {
    #                 velocityField:'vectorField',
    #                 pressureField:'scalarField'
    #             }
    #             self._parameters = {
    #                 dynViscosity:'scalar'
    #             }
    #         def show(self):
    #             print("IncompressibleFlow:\n\tvariables   {o._variables}\n\tparameters  {o._parameters}\n".format(o=self))
    #
    # class Mesh:
    #     def __init__(self, geometry, kernel):
    #         self._geometry = geometry
    #         self._globalrefinement = 1

    # this class must contain a list of available flow models and boundaries. maybe ask kernel for it
    class ContinuumProperties:
        def __init__(self, flowmodels, geometry):
            # these are too many quite similar dicts
            self._variableFields = {}
            self._parameterFields = {}
            self._boundaries = {}
            self._variableBoundaries = {}

            self.collocateFields(flowmodels)
            self.defineBoundariesOnFields(geometry)

        def collocateFields(self, flowmodels):
            for fm in flowmodels:
                for v,t in fm._variables.items():
                    self._variableFields[v] = t
                # self._variableFields += fm._variables
                for p,t in fm._parameters.items():
                    self._parameterFields[p] =t
            #
            # # removing dublicate field entries
            # self._variableFields = list(set(self._variableFields))
            # self._parameterFields = list(set(self._parameterFields))

            # parameters that also appear in variable list are removed
            for v in self._variableFields:
                self._parameterFields.pop(v, None)
            # for p in self._parameterFields:
            #     if p in self._variableFields:
            #         self._parameterFields.remove(p)

        def defineBoundariesOnFields(self, geometry):
            for b in geometry._boundaries:
                self._boundaries[b] = 'vNeumann'        #default pointer, this should be a pointer to a bc object

            for v in self._variableFields:
                self._variableBoundaries[v] = self._boundaries.copy()

        def setBoundaryCondition(self, field, boundaryName, boundaryType, kwargs):
            self._variableBoundaries[field][boundaryName] = boundaryType

        def show(self):
            print("variable fields ", self._variableFields)
            print("parameter fields ", self._parameterFields)
            for v,b in self._variableBoundaries.items():
                print(v)
                print(b)
            print("\n")
