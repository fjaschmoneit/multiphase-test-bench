import Fields, LinearEquationSystems
import PrimitiveFields
import FVM

class ScalarConvectionDiffusion():
    def __init__(self, scalarFieldName='s', velocityFieldName='U', diffusionCoefficientName='D', sourceFieldName='R'):
        self._variables = {
            scalarFieldName: 'scalarField'
        }
        self._parameters = {
            velocityFieldName: 'vectorField',
            sourceFieldName : 'scalarField',
            diffusionCoefficientName: 'scalar'
        }
        self._properties = {
            'diffusionCoefficient': None,
            'sourceFieldName': sourceFieldName,
            'scalarFieldName': scalarFieldName,
            'velocityFieldName': velocityFieldName
        }
        self.availableBoundaryConditions = None

    def showContinuumProperties(self):
        print("ScalarConvectionDiffusion:\n\t", self._properties)

    def updateLinearEquationSystems(self, mesh, fields):

        # select correct field:
        depField = fields[self._properties['scalarFieldName']]
        sourceField = fields[self._properties['sourceFieldName']]
        velocityField = fields[self._properties['velocityFieldName']]

        # definings convective flux, vectorField
        F = velocityField

        # defining diffusive flux, vectorField
        A = mesh.getFaceAreas()
        D = mesh.getInverseCellDistances() * self._properties['diffusionCoefficient'] * A

        # defining source, scalarField
        S = sourceField * mesh._uniformSpacing

        depField._A, depField._b = LinearEquationSystems.createCoefficientMatrix(mesh,F,D,S, depField)

class IncompressibleFlow():
    def __init__(self, velocityField='U', pressureField='p', dynViscosity='nu'):
        self._variables = {
            velocityField: 'vectorField',
            pressureField: 'scalarField'
        }
        self._parameters = {
            dynViscosity: 'scalar'
        }
        self._properties = {
            'dynamic viscosity [m^2/s]': 1e-6,
            'pressureFieldName': 'p',
            'velocityFieldName': 'U',
        }

    # creating staggered fields for u,v and
    def createLinearEquationSystems(self, mesh, fields):
        # select correct field:
        velocityField = fields[self._properties['velocityFieldName']]
        sourceField = fields[self._properties['sourceFieldName']]
        pressureField = fields[self._properties['velocityFieldName']]

        fluxes = self.createMatrixCoefficients(
            mesh=mesh,
            phi=depField,
            U=velocityField,
            diffCoeff=self._properties['diffusionCoefficient'],
            sourceField=sourceField)

        depField._A, depField._b = LinearEquationSystems.createCoefficientMatrix(*fluxes)

    def createMatrixCoefficients(self, mesh, phi, diffCoeff, U=None, rho=None, sourceField=None):
        ### returning coefficient vectors for subsequent coefficient matrix assembly and coeff vector

        #returns either zero, if we have pure diffusion, otherwise a parameterFaceField
        F = self.createConvectiveFluxes( rho,U  )
        D = mesh.getInverseCellDistances() * diffCoeff
        A = mesh.getFaceAreas()

        # directional coefficient matrices are cellFields, very inefficient
        a_e = -Fields.parameterCellField(mesh=mesh, primitiveField=(D*A).e - 0.5*F.e )
        a_w = -Fields.parameterCellField(mesh=mesh, primitiveField=(D*A).w + 0.5*F.w )
        a_n = -Fields.parameterCellField(mesh=mesh, primitiveField=(D*A).n - 0.5*F.n )
        a_s = -Fields.parameterCellField(mesh=mesh, primitiveField=(D*A).s + 0.5*F.s )

        s_u = sourceField * mesh._uniformSpacing

        a_p = -(a_e + a_w + a_n + a_s)

        # fixing boundary conditions
        if phi._boundary['top'] == 'zeroGradient':
            a_p.bn += a_n.bn
            a_n.bn = 0
        else:
            Tn = phi._boundary['top']
            s_u.bn -= a_n.bn * Tn
            a_n.bn = 0

        if phi._boundary['bottom'] == 'zeroGradient':
            a_p.bs += a_s.bs
            a_s.bs = 0
        else:
            Ts = phi._boundary['bottom']
            s_u.bs -= a_s.bs * Ts
            a_s.bs = 0

        if phi._boundary['left'] == 'zeroGradient':
            a_p.bw += a_w.bw
            a_w.bw = 0
        else:
            Tw = phi._boundary['left']
            s_u.bw -= a_w.bw * Tw
            a_w.bw = 0

        if phi._boundary['right'] == 'zeroGradient':
            a_p.be += a_e.be
            a_e.be = 0
        else:
            Te = phi._boundary['right']
            s_u.be -= a_e.be * Te
            a_e.be = 0

        return a_e._raw, a_w._raw, a_n._raw, a_s._raw, a_p._raw, s_u._raw


    def show(self):
        print("IncompressibleFlow:\n\tvariables   {o._variables}\n\tparameters  {o._parameters}\n".format(o=self))


