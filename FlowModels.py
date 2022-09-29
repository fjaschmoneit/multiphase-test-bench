import Interpolation

class ScalarConvectionDiffusion():

    def __init__(self, depVariableName, **fieldnames):
        self._mesh = None
        self._depVariableName = depVariableName
        self._velocityFieldName = fieldnames['velocityFieldName']
        self._diffusionCoeffName = fieldnames['diffusionCoeffName']
        self._sourceField_p = None
        self._sourceField_c = None
        self._convFluxes = None
        self._diffFluxes = None

    # I should have global instances of these fields
    def initialize(self, mesh, fieldReg, fieldFlowModelLink):
        self._mesh = mesh
        fGov = fieldReg['governor']

        self._sourceField_p = fGov.newField(type='scalarCV')
        self._sourceField_c = fGov.newField(type='scalarCV')

        #self._sourceField = fGov.sourceField(shape=fGov.shapeCVField)
        self._convFluxes = fGov.newFaceField(shape_u=fGov.shapeFaces_u, shape_v=fGov.shapeFaces_v)
        self._diffFluxes = fGov.newFaceField(shape_u=fGov.shapeFaces_u, shape_v=fGov.shapeFaces_v)

        self.initializeFields(fieldReg)
        self.linkDepFieldToModel(fieldFlowModelLink)

    # def updateSourceField(self, fieldReg):
    #     constHeatSource = fieldReg[self._depVariableName]._source
    #     self._sourceField_c.fill( constHeatSource )
    #     self._sourceField_p.fill( constHeatSource )
    #     #self._sourceField.setConstantSource( constHeatSource*self._mesh._uniformSpacing )

    def linkDepFieldToModel(self, fieldFlowModelLink):
        fieldFlowModelLink[self._depVariableName] = self

    def initializeFields(self, fieldReg):
        # adding neccessary field to registry
        # check if fields are not already defined

        fGov = fieldReg['governor']
        fieldReg[self._depVariableName] = fGov.newField(type='scalarCV', includeGhostNodes=True)
        #fieldReg[self._depVariableName] = fGov.newScalarField(shape=fGov.shapeVarScalarField)
        fieldReg[self._velocityFieldName] = fGov.newFaceField(shape_u=fGov.shapeFaces_u, shape_v=fGov.shapeFaces_v)
        fieldReg[self._diffusionCoeffName] = 0.0

    def setConstSource(self, value):
        self._sourceField_c.fill(value*self._mesh._uniformSpacing)

    def setDerichlet(self, loc, value):
        D = self._diffFluxes
        F = self._convFluxes
        Sc = self._sourceField_c
        Sp = self._sourceField_p

        if loc == 'left':
            Sc.bw += ( 2.0 * D.u.bw + F.u.bw )* value
            Sp.bw += 2.0 * D.u.bw + F.u.bw
            D.u.bw = 0.0
            F.u.bw = 0.0
        elif loc == 'right':
            Sc.be += ( 2.0 * D.u.be - F.u.be) * value
            Sp.be += 2.0 * D.u.be - F.u.be
            D.u.be = 0.0
            F.u.be = 0.0
        elif loc == 'top':
            Sc.bn += ( 2.0 * D.v.bn - F.v.bn) * value
            Sp.bn +=  2.0 * D.v.bn - F.v.bn
            D.v.bn = 0.0
            F.v.bn = 0.0
        elif loc == 'bottom':
            Sc.bs += ( 2.0 * D.v.bs + F.v.bs) * value
            Sp.bs += 2.0 * D.u.bs + F.v.bs
            D.v.bs = 0.0
            F.v.bs = 0.0

    def setVonNeumann(self, loc):
        D = self._diffFluxes
        # Sc = self._sourceField.Sc
        # Sp = self._sourceField.Sp

        if loc == 'left':
            D.u.bw = 0.0
        elif loc == 'right':
            D.u.be = 0.0
        elif loc == 'top':
            D.v.bn = 0.0
        elif loc == 'bottom':
            D.v.bs = 0.0

    def correctBCs(self, fieldReg):
        for loc, type in fieldReg[self._depVariableName]._boundary.items():
            if type == 'zeroGradient':
                self.setVonNeumann(loc)
            else:
                value = type  # I need a better interface
                self.setDerichlet(loc, value)

    # this is the most important FM method. The others could possibly be base class methods
    def updateFluxes(self, fieldReg):
        self._convFluxes = fieldReg[self._velocityFieldName] * self._mesh.calcFaceAreas(fieldReg['governor'])
        self._diffFluxes = fieldReg['invCellDist'] * fieldReg[self._diffusionCoeffName] * self._mesh.calcFaceAreas(fieldReg['governor'])

# make inherit from some base flowmodel class
class IncompressibleMomentumComp:

    def __init__(self, depVariableName, **fieldnames):
        self._mesh = None
        self._sourceField = a
        self._depVariableName = depVariableName
        self._orientation = fieldnames['orientation']
        self._otherVelFieldName = fieldnames['otherVelocityFieldName']
        self._pressureFieldName = fieldnames['pressureFieldName']
        self._kinViscosity = fieldnames['kinViscosityName']

    # I should have global instances of these fields
    def initialize(self, mesh, fieldReg, fieldFlowModelLink):
        self._mesh = mesh
        fGov = fieldReg['governor']
        self._convFluxes = fGov.newVectorField(shape_u=fGov.shapeCVField, shape_v=fGov.shapeVerticesInternal_u)
        self._diffFluxes = fGov.newVectorField(shape_u=fGov.shapeCVField, shape_v=fGov.shapeVerticesInternal_u)
        self._sourceField = fGov.sourceField(shape=fGov.shapeFacesInternal_u)

        self.initializeFields(fieldReg)
        self.linkDepFieldToModel(fieldFlowModelLink)

    # do I need a source field?
    def updateSourceField(self, fieldReg):
        self._sourceField.Sc.fill(0.0)
        self._sourceField.Sp.fill(0.0)

    # this could be a base class
    def linkDepFieldToModel(self, fieldFlowModelLink):
        fieldFlowModelLink[self._depVariableName] = self

    def initializeFields(self, fieldReg):
        # adding neccessary field to registry
        # check if fields are not already defined
        fGov = fieldReg['governor']

        fieldReg[self._depVariableName] = self.linkOrCreateField(self._depVariableName, fieldReg, fGov.shapeFaces_u)
        fieldReg[self._otherVelFieldName] = self.linkOrCreateField(self._otherVelFieldName, fieldReg, fGov.shapeFaces_v)
        fieldReg[self._pressureFieldName] = self.linkOrCreateField(self._pressureFieldName, fieldReg, fGov.shapeCVField)
        fieldReg[self._kinViscosity] = 0.0



    def correctBCs(self, fieldReg):
        for loc, type in fieldReg[self._depVariableName]._boundary.items():
            if type == 'zeroGradient':
                self.setVonNeumann(loc)
            else:
                value = type  # I need a better interface
                self.setDerichlet(loc, value)

    # this is the most important FM method. The others could possibly be base class methods
    def updateFluxes(self, fieldReg):
        u = fieldReg[self._depVariableName]
        v = fieldReg[self._otherVelFieldName]
        rCellDist_u = fieldReg['invCellDist'].u
        rCellDist_v = fieldReg['invCellDist'].v

        self._convFluxes.u = 0.5*( u.east + u.west )
        self._convFluxes.v = 0.5 *(v.east + v.west)

        self._diffFluxes.u = 0.5*( rCellDist_u.east + rCellDist_u.west ) * fieldReg[self._kinViscosity]
        self._diffFluxes.v = 0.5*( rCellDist_v.east + rCellDist_v.west )* fieldReg[self._kinViscosity]

    # this could be a base method
    def linkOrCreateField(self, name, fieldReg, shape):
        if name not in fieldReg:
            fGov = fieldReg['governor']
            return fGov.newScalarField(shape=shape)
        else:
            return fieldReg[name]

    def setDerichlet(self, loc, value):
        D = self._diffFluxes
        F = self._convFluxes
        Sc = self._sourceField.Sc
        Sp = self._sourceField.Sp

        if loc == 'left':
            Sc.bw += ( 2.0 * D.u.bw + F.u.bw )* value
            Sp.bw += 2.0 * D.u.bw + F.u.bw
            D.u.bw = 0.0
            F.u.bw = 0.0
        elif loc == 'right':
            Sc.be += ( 2.0 * D.u.be - F.u.be) * value
            Sp.be += 2.0 * D.u.be - F.u.be
            D.u.be = 0.0
            F.u.be = 0.0
        elif loc == 'top':
            Sc.bn += ( 2.0 * D.v.bn - F.v.bn) * value
            Sp.bn +=  2.0 * D.v.bn - F.v.bn
            D.v.bn = 0.0
            F.v.bn = 0.0
        elif loc == 'bottom':
            Sc.bs += ( 2.0 * D.v.bs + F.v.bs) * value
            Sp.bs += 2.0 * D.u.bs + F.v.bs
            D.v.bs = 0.0
            F.v.bs = 0.0

    def setVonNeumann(self, loc):
        D = self._diffFluxes
        # Sc = self._sourceField.Sc
        # Sp = self._sourceField.Sp

        if loc == 'left':
            D.u.bw = 0.0
        elif loc == 'right':
            D.u.be = 0.0
        elif loc == 'top':
            D.v.bn = 0.0
        elif loc == 'bottom':
            D.v.bs = 0.0









