import FlowModels, Mesh, Fields, Differentiation, LinearEquationSystems

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

        #self._continuumProperties = None
        #self._discreteProperties = None

        self._variableFields = {}
        self._coefficientFields = {}
        self._scalarCoefficients = {}
        #self._linearSystems = {}  # in fields

        self.createfields()   # fields should be created in compile()

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

    def deleteFields(self):
        self._variableFields.clear()
        self._coefficientFields.clear()

    def getFields(self):
        return self._variableFields

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


    def compile(self):
        # creating fields and their corresponding matrix equations

        #self.createfields() are already created. change this

        # creates the linear equation systems
        for flowmodel in self._flowmodels:
            flowmodel.createLinearEquationSystems(mesh=self._mesh, fields=self._variableFields)

        self._isCompiled = True

    # also defining the solution algorithm, i.e. sequence of SOLVING variable fields and UPDATING coefficient fields
    # def initialize(self, flowmodels, continuumProperties, mesh):
    #     self._mesh = mesh
    #     self._flowmodels = flowmodels
    #     self._continuumProperties = continuumProperties
    #
    #     self.createFields()  #are these fields not only temporary, since I calculate in matrices all the time?
    #     #self.createLinearEquationSystems()
    #     # self.defineIterationSquence


    # this class calls backend methods
# class simulation():
#
#     def __init__(self, parameters):
#         self._mesh = None
#         self._variableFields = {}
#         self._coefficientFields = {}
#         self._linearSystems = {}
#
#         self.createMesh(parameters['mesh'])
#         self.createFields(parameters['fields'])  #are these fields not only temporary, since I calculate in matrices all the time?
#         self.createLinearEquationSystems(parameters['flowModels'])


        # creating a linear equation system for every field variable:
        # how to include bcs here?

            # fluxes = flowmodel.ScalarConvectionDiffusion.createMatrixCorefficients(mesh=self._mesh, field=self._variableFields, diffCoeff=10)
            # field._A, field._b = LinearEquationSystems.createCoefficientMatrix(*fluxes)

        # for flowmodel in self._flowmodels:
        #     FlowModels.createLinearEquationSystems(flowmodel, self._mesh, self._variableFields)


            # if isinstance(fm, FlowModels.ScalarConvectionDiffusion):
            #     fluxes = FlowModels.scalarDiffusion(mesh=self._mesh, field=field, diffCoeff=10)
            #     A, b = LinearEquationSystems.createCoefficientMatrix(*fluxes)
            #     self._linearSystems[field] = {
            #         'A' : A,
            #         'b' : b
            #     }

        # for field in self._variableFields:
        #     if flowModelPrms[field]['type'] == 'scalarDiffusion':
        #         fluxes = FlowModels.scalarDiffusion(mesh=self._mesh, field=field, diffCoeff=10)
        #         A, b = LinearEquationSystems.createCoefficientMatrix(*fluxes)
        #         self._linearSystems[field] = {
        #             'A' : A,
        #             'b' : b
        #         }
        #     else:
        #         print("error: no flow model chosen for field ", field)

        # for fieldName in fieldPrmsDict:
        #     if fieldPrmsDict[fieldName]['type'] == 'scalarField':
        #         if fieldPrmsDict[fieldName]['modifier'] == 'iterative':
        #             self._variableFields[fieldName] = Fields.parameterCellField( mesh=self._mesh,
        #                                                                 value=fieldPrmsDict[fieldName]['initialValue'] )
        #         elif fieldPrmsDict[fieldName]['modifier'] == 'algebraic':
        #             self._coefficientFields[fieldName] = Fields.parameterCellField( mesh=self._mesh,
        #                                                                 value=fieldPrmsDict[fieldName]['initialValue'] )
        #         else:
        #             print("Error: modifier in field not recognized")
        #     elif fieldPrmsDict[fieldName]['type'] == 'constant':
        #         self._coefficientFields[fieldName] = fieldPrmsDict[fieldName]['initialValue']
        #     else:
        #         print("unknown field type")

    # #should not be a kernel method
    # def createMesh(self, meshPrms):
    #     if meshPrms['type'] == '2DCartesian':
    #         self._mesh = Mesh.cartesian2D(
    #             len_x=meshPrms['length_x'],
    #             len_y=meshPrms['length_y'],
    #             res=meshPrms['resolution'] )
    #     else:
    #         print("incorrect mesh type")

    def execute(self):
        if self._isCompiled:
            print("executing")
            return 1
        else:
            print("simulation object not compiled")


    def display(self, field):
        Fields.drawCellField(field)


    def prepSim(self, simPrms):
        self.createMesh(simPrms['mesh'])
        #self._mesh = Mesh.cartesian2D( len_x=0.5, len_y=0.3, res=100 )       #doesnt work for higher res or len_y
        #print("reciprocal mesh dist\n", self.mesh.invCellDist.ff_y)

        # parameter field without boundary conditions
        self._fields['gamma'] = Fields.parameterCellField( mesh=self._mesh, value=0 )
        self._fields['gamma'].fillWithConsecutiveValues()
        #print("gamma cellField:\n",gamma.internal)

        # this should become a constructor in a class
        # a parameter Face Field
        #self.fieldDict['gamma_f'] = Interpolation.centralDifferencing(self.fieldDict['gamma'], self.mesh)
        #
        # print("centralDiff faceField entire x\n",gamma_f.ff_x)
        # print("centralDiff faceField entire y\n",gamma_f.ff_y)
        self._T = Fields.parameterCellField( mesh=self._mesh, value=0 )

        self._T.fillWithRandomIntegers()

        #print("T = \n", T.raw)

        gradT = Differentiation.grad(cellField=self._T, mesh=self._mesh)
        # print("grad T _x\n", gradT.ff_x)
        # print("grad T _y\n", gradT.ff_y)

        fluxes = FlowModels.scalarDiffusion(mesh=self._mesh, field=self._T, diffCoeff=10)

        self._A, self._b = LinearEquationSystems.createCoefficientMatrix(*fluxes)

#
# # I basically need two kinds of fields, a coefficient field, and a quantity field
#
# print("cell field: \n", T.field)   # this field is evaluated, therefore I also need to specify boundary conditions
# print("gradient inner: \n", gradT[0].internal)
# print("gradient boundary: \n", gradT[0].boundary)   # relates to boundary condition
# print("gradient entire: \n", gradT[0].field)



#T.setTransportEuqation(name='diffusion', gamma='10')
# T.setBoundary(name='left', type='fixedValue', value=100 )
# T.setBoundary(name='right', type='fixedValue', value=500 )
# T.createCoefficientMatrix()
# T.solve()



# exapmle from book, just for reference
# T = np.arange(5)
# coeffMatrixObj = bookExample.coeffMatrix(dim = len(T))
# coeffMatrixObj.setCoeffsMatrix()
# coeffMatrixObj.setCoeffVector(100,500)
# A = coeffMatrixObj.A
# b = coeffMatrixObj.b
# T = np.linalg.solve(A, b)
# print("book solution")
# print(T)

