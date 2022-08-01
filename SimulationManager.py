import FlowModels, Mesh, Fields, Differentiation, LinearEquationSystems

class simulation():

    def __init__(self):
        self.fieldDict = {}

    def prepSim(self):
        self.mesh = Mesh.cartesian2D( len_x=0.5, len_y=0.3, res=100 )       #doesnt work for higher res or len_y
        #print("reciprocal mesh dist\n", self.mesh.invCellDist.ff_y)

        # parameter field without boundary conditions
        self.fieldDict['gamma'] = Fields.parameterCellField( mesh=self.mesh, value=0 )
        self.fieldDict['gamma'].fillWithConsecutiveValues()
        #print("gamma cellField:\n",gamma.internal)

        # this should become a constructor in a class
        # a parameter Face Field
        #self.fieldDict['gamma_f'] = Interpolation.centralDifferencing(self.fieldDict['gamma'], self.mesh)
        #
        # print("centralDiff faceField entire x\n",gamma_f.ff_x)
        # print("centralDiff faceField entire y\n",gamma_f.ff_y)
        self.T = Fields.parameterCellField( mesh=self.mesh, value=0 )

        self.T.fillWithRandomIntegers()

        #print("T = \n", T.raw)

        gradT = Differentiation.grad(cellField=self.T, mesh=self.mesh)
        # print("grad T _x\n", gradT.ff_x)
        # print("grad T _y\n", gradT.ff_y)

        fluxes = FlowModels.scalarDiffusion(mesh=self.mesh, field=self.T, diffCoeff=10)

        self.A, self.b = LinearEquationSystems.createCoefficientMatrix(*fluxes)


    def execute(self):
        x = LinearEquationSystems.solveLinearSystem(self.A,self.b)
        self.T.raw[:,:] = x.reshape(self.T.ny, self.T.nx)
        Fields.drawCellField(self.T)



#
# A = np.ones((3,2))
# print(A)
# B = np.reshape(A, 6)
# print(B)
#
# A[1][1] = 5
# print(A)
# print(B)




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

