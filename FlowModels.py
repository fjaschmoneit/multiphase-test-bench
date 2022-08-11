import Fields, LinearEquationSystems

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
            'diffusionCoefficient': 1,
            'sourceFieldName': sourceFieldName,
            'scalarFieldName': scalarFieldName,
            'velocityFieldName': velocityFieldName
        }
        self.availableBoundaryConditions = None

    def showContinuumProperties(self):
        print("ScalarConvectionDiffusion:\n\t", self._properties)

    def createLinearEquationSystems(self, mesh, fields):

        # select correct field:
        depField = fields[self._properties['scalarFieldName']]
        sourceField = fields[self._properties['sourceFieldName']]

        fluxes = self.createMatrixCorefficients(
            mesh=mesh,
            field=depField,
            diffCoeff=self._properties['diffusionCoefficient'],
            sourceField=sourceField)

        depField._A, depField._b = LinearEquationSystems.createCoefficientMatrix(*fluxes)

    def createMatrixCorefficients(self, mesh, field, diffCoeff, sourceField):
        ### returning coefficient vectors for subsequent coefficient matrix assembly and coeff vector

        faceFluxes = Fields.parameterFaceField(mesh=mesh)
        gamma = Fields.parameterFaceField(mesh=mesh, value=diffCoeff)

        # a mesh method
        faceArea = Fields.parameterFaceField(mesh=mesh, value=1.0)
        #cellVolumes = mesh.getCellVolumes()

        fax = faceArea.entries_EW
        fay = faceArea.entries_NS
        gammax = gamma.entries_EW
        gammay = gamma.entries_NS
        qx = mesh.invCellDist.entries_EW
        qy = mesh.invCellDist.entries_NS

        # defining internal values
        faceFluxes.entries_EW = fax * gammax * qx
        faceFluxes.entries_NS = fay * gammay * qy

        # directional coefficient matrices are cellFields, very inefficient
        a_e = -Fields.parameterCellField(mesh=mesh, primitiveField=faceFluxes.e)
        a_w = -Fields.parameterCellField(mesh=mesh, primitiveField=faceFluxes.w)
        a_n = -Fields.parameterCellField(mesh=mesh, primitiveField=faceFluxes.n)
        a_s = -Fields.parameterCellField(mesh=mesh, primitiveField=faceFluxes.s)

        s_u = sourceField * mesh._uniformSpacing *1.0

        a_p = -(a_e + a_w + a_n + a_s)

        # fixing boundary conditions
        if field._boundary['top'] == 'zeroGradient':
            a_p.bn += a_n.bn
            a_n.bn = 0
        else:
            Tn = field._boundary['top']
            s_u.bn -= a_n.bn * Tn
            a_n.bn = 0

        if field._boundary['bottom'] == 'zeroGradient':
            a_p.bs += a_s.bs
            a_s.bs = 0
        else:
            Ts = field._boundary['bottom']
            s_u.bs -= a_s.bs * Ts
            a_s.bs = 0

        if field._boundary['left'] == 'zeroGradient':
            a_p.bw += a_w.bw
            a_w.bw = 0
        else:
            Tw = field._boundary['left']
            s_u.bw -= a_w.bw * Tw
            a_w.bw = 0

        if field._boundary['right'] == 'zeroGradient':
            a_p.be += a_e.be
            a_e.be = 0
        else:
            Te = field._boundary['right']
            s_u.be -= a_e.be * Te
            a_e.be = 0


        return a_e._raw, a_w._raw, a_n._raw, a_s._raw, a_p._raw, s_u._raw


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

    def show(self):
        print("IncompressibleFlow:\n\tvariables   {o._variables}\n\tparameters  {o._parameters}\n".format(o=self))


