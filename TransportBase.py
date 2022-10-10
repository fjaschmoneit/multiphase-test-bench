import Fields

class transportBase:

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem, depFieldShape):
        self._mesh = mesh
        self._linSystem = linSystem
        self._fieldReg = fieldReg
        self._fc = fieldCreator
        self._depFieldShape = depFieldShape
        self._depField = self._fc.newField(shape=self._depFieldShape)
        self._diffusionCoefficient = None
        self._convFluxes = None
        self._diffFluxes = None
        self._sourceField_c = None
        self._sourceField_p = None

        self.createConfDiffFluxes()

    def createConfDiffFluxes(self):
        (shape_u, shape_v) = self._fc.getFluxShapesFromCVShape(self._depFieldShape)
        self._convFluxes = Fields.fieldContainer(
            u=self._fc.newField(shape=shape_u),
            v=self._fc.newField(shape=shape_v)
        )

        self._diffFluxes = Fields.fieldContainer(
            u=self._fc.newField(shape=shape_u),
            v=self._fc.newField(shape=shape_v)
        )

    def setSourceField(self, value):
        self._outerSource = self._fc.newField(shape=self._depFieldShape, value=value)

    def setDiffusionCoefficient(self, value):
        self._diffusionCoefficient = value

    def setConvectionField(self, **kwargs):
        (const_u, const_v) = kwargs.get('values', (None, None))

        (shape_u, shape_v) = self._fc.getFluxShapesFromCVShape(self._depFieldShape)

        self._convectionField = kwargs.get('field', Fields.fieldContainer(
            u=self._fc.newField(shape=shape_u, value=const_u),
            v=self._fc.newField(shape=shape_v, value=const_v)
        ))

    def getDepField(self):
        return self._depField

    def correctBCs(self):
        self._sourceField_p = self._fc.newField(shape=self._depFieldShape, value=0.0)

        # self._sourceField_c = self._fc.newField(data=self._outerSource.data)  # does not work??
        self._sourceField_c = self._fc.newField(shape=self._depFieldShape, value=0.0)
        for loc, typeValueTupel in self._depField._boundary.items():
            (type, value) = typeValueTupel
            if type == 'zeroGradient':
                self.setVonNeumann(loc)
            else:
                self.setDerichlet(loc, value)

    def updateLinSystem(self):
        self._linSystem.update(F=self._convFluxes, D=self._diffFluxes, Sc=self._sourceField_c, Sp=self._sourceField_p)

    def solve(self):
        self._depField.data = self._linSystem.solve()

    def setVonNeumann(self, loaction):
        pass

    def setDerichlet(self, location, value):
        pass


# make inherit from some base flowmodel class
class StaggeredTransport:

    def __init__(self, depVariableName, **fieldnames):
        self._depVariableName = depVariableName
        self._otherVelFieldName = fieldnames['otherVelocityFieldName']
        self._pressureFieldName = fieldnames['pressureFieldName']
        self._kinViscosity = fieldnames['kinViscosityName']
        self._boundaryDict = {}

    def initializeFlowModel(self, mesh, fieldReg, fieldFlowModelLink):
        self._mesh = mesh
        fGov = fieldReg['governor']

        # these container dimensions are only valid for momentum U component
        # must make abstraction
        self._convFluxes = Fields.fieldContainer(
            u=fGov.newField(type='staggered_u', value=0.0),
            v=fGov.newField(type='vertex', value=0.0)
        )

        self._diffFluxes = Fields.fieldContainer(
            u=fGov.newField(type='staggered_u', value=0.0),
            v=fGov.newField(type='vertex', value=0.0)
        )

        self._sourceField_c = fGov.newField(type='faces_u')
        self._sourceField_p = fGov.newField(type='faces_u')

        self.initializeFields(fieldReg)
        self.linkDepFieldToModel(fieldFlowModelLink)

    # this could be a base class
    def linkDepFieldToModel(self, fieldFlowModelLink):
        fieldFlowModelLink[self._depVariableName] = self

    def initializeFields(self, fieldReg):
        # adding neccessary field to registry
        # check if fields are not already defined

        fieldReg[self._depVariableName] = self.linkOrCreateField(self._depVariableName, fieldReg, type='faces_u')
        fieldReg[self._otherVelFieldName] = self.linkOrCreateField(self._otherVelFieldName, fieldReg, type='faces_v')
        fieldReg[self._pressureFieldName] = self.linkOrCreateField(self._pressureFieldName, fieldReg, type='scalarCV')
        fieldReg[self._kinViscosity] = 0.0

    def correctBCs(self, fieldReg):
        for loc, typeValueTupel in fieldReg[self._depVariableName]._boundary.items():
            (type, value) = typeValueTupel
            if type == 'zeroGradient':
                self.setVonNeumann(loc)
            else:
                self.setDerichlet(loc, value)

    # this is the most important FM method. The others could possibly be base class methods
    def updateFluxes(self, fieldReg):

        u = fieldReg[self._depVariableName]
        v = fieldReg[self._otherVelFieldName]
        rCellDist_u = fieldReg['invCellDist'].u
        rCellDist_v = fieldReg['invCellDist'].v

        self._convFluxes.u.internal_u = 0.5 * (u.east + u.west)
        self._convFluxes.v.internal_u = 0.5 * (v.east + v.west)

        self._diffFluxes.u.internal_u = 0.5*( rCellDist_u.east + rCellDist_u.west ) * fieldReg[self._kinViscosity]
        self._diffFluxes.v.internal_u = 0.5*( rCellDist_v.east + rCellDist_v.west )* fieldReg[self._kinViscosity]

    def updateSourceField(self, fieldReg, mesh):
        p = fieldReg[self._pressureFieldName]
        faceAreas_u = mesh.calcFaceAreas(fieldReg['governor']).u
        gradP_u = Differentiation.grad_u(p, fieldReg)
        self._sourceField_c.internal_u = -gradP_u*faceAreas_u.internal_u

    # this could be a base method
    def linkOrCreateField(self, name, fieldReg, type):
        if name not in fieldReg:
            fGov = fieldReg['governor']
            return fGov.newField(type=type)
        else:
            return fieldReg[name]

    def setDerichlet(self, loc, value):
        D = self._diffFluxes
        F = self._convFluxes
        Sc = self._sourceField_c
        Sp = self._sourceField_p

        if loc == 'left':
            Sc.bw += 1e15*value
            Sp.bw += 1e15
        elif loc == 'right':
            Sc.be += 1e15 * value
            Sp.be += 1e15
        elif loc == 'top':
            Sc.bn += 1e15 * value
            Sp.bn += 1e15
        elif loc == 'bottom':
            Sc.bs += 1e15 * value
            Sp.bs += 1e15

    def setVonNeumann(self, loc):
        D = self._diffFluxes
        F = self._convFluxes
        Sc = self._sourceField_c
        Sp = self._sourceField_p

        #I just set these fluxes to zero. But I should set them to their intenal neighbors value
        if loc == 'left':
            D.u.bw = 0.0
            F.u.bw = 0.0
        elif loc == 'right':
            D.u.be = 0.0
            F.u.be = 0.0
        elif loc == 'top':
            D.v.bn = 0.0
            F.v.bn = 0.0
        elif loc == 'bottom':
            D.v.bs = 0.0
            F.v.bs = 0.0

    # implement difference schemes here
    #def getTotalFluxes(self):

    def getTotalFluxes_e(self):
        return -(self._diffFluxes.u.east - 0.5*self._convFluxes.u.east)

    def getTotalFluxes_w(self):
        return -(self._diffFluxes.u.west + 0.5 * self._convFluxes.u.west)

    def getTotalFluxes_n(self):
        return -(self._diffFluxes.v.north - 0.5 * self._convFluxes.v.north)

    def getTotalFluxes_s(self):
        return -(self._diffFluxes.v.south - 0.5 * self._convFluxes.v.south)
