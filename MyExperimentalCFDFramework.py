import FlowModels, Mesh, Fields, Differentiation, LinearEquationSystems

# this class calls backend methods
class simulation():

    def __init__(self, parameters):
        self._mesh = None
        self._variableFields = {}
        self._coefficientFields = {}
        self._linearSystems = {}

        self.createMesh(parameters['mesh'])
        self.createFields(parameters['fields'])  #are these fields not only temporary, since I calculate in matrices all the time?
        self.createLinearEquationSystems()

    def createLinearEquationSystems(self):
        for field in self._variableFields:
            fluxes = FlowModels.scalarDiffusion(mesh=self._mesh, field=field, diffCoeff=10)
            A, b = LinearEquationSystems.createCoefficientMatrix(*fluxes)
            self._linearSystems[field] = {
                'A' : A,
                'b' : b
            }

    def createFields(self, fieldPrmsDict):
        for fieldName in fieldPrmsDict:
            if fieldPrmsDict[fieldName]['type'] == 'scalarField':
                if fieldPrmsDict[fieldName]['modifier'] == 'iterative':
                    self._variableFields[fieldName] = Fields.parameterCellField( mesh=self._mesh,
                                                                        value=fieldPrmsDict[fieldName]['initialValue'] )
                elif fieldPrmsDict[fieldName]['modifier'] == 'algebraic':
                    self._coefficientFields[fieldName] = Fields.parameterCellField( mesh=self._mesh,
                                                                        value=fieldPrmsDict[fieldName]['initialValue'] )
                else:
                    print("Error: modifier in field not recognized")
            elif fieldPrmsDict[fieldName]['type'] == 'constant':
                self._coefficientFields[fieldName] = fieldPrmsDict[fieldName]['initialValue']
            else:
                print("unknown field type")

    def createMesh(self, meshPrms):
        if meshPrms['type'] == '2DCartesian':
            self._mesh = Mesh.cartesian2D(
                len_x=meshPrms['length_x'],
                len_y=meshPrms['length_y'],
                res=meshPrms['resolution'] )
        else:
            print("incorrect mesh type")


    def solve(self, field):
        x = LinearEquationSystems.solveLinearSystem(self._linearSystems[field]['A'], self._linearSystems[field]['b'])
        self._variableFields[field].raw[:, :] = x.reshape(self._variableFields[field].ny, self._variableFields[field].nx)
        #self._T.raw[:, :] = x.reshape(self._T.ny, self._T.nx)

    def display(self, fieldName):
        Fields.drawCellField(self._variableFields[fieldName])


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

