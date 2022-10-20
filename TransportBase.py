import Fields

class transportBase:

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem, depFieldShape):
        self._mesh = mesh
        self._linSystem = linSystem
        self._fieldReg = fieldReg
        self._fc = fieldCreator
        self._depFieldShape = depFieldShape
        self._diffusionCoefficient = None
        self._convFluxes = None
        self._diffFluxes = None

        self.depField = self._fc.newField(shape=self._depFieldShape, value=0.0)  # why is this not a simple data field?
        self._sourceField_c = Fields.newDataField(shape=self._depFieldShape, value=0.0)
        self.a_e = Fields.newDataField(shape=self._depFieldShape, value=0.0)
        self.a_w = Fields.newDataField(shape=self._depFieldShape, value=0.0)
        self.a_s = Fields.newDataField(shape=self._depFieldShape, value=0.0)
        self.a_n = Fields.newDataField(shape=self._depFieldShape, value=0.0)
        self.a_p = Fields.newDataField(shape=self._depFieldShape, value=0.0)

    def reset(self):
#        self.depField.data.fill(0.0)
        self._sourceField_c.fill(0.0)
        self.a_e.fill(0.0)
        self.a_w.fill(0.0)
        self.a_s.fill(0.0)
        self.a_n.fill(0.0)
        self.a_p.fill(0.0)

    def setSourceField(self, value):
        self._sourceField_c.fill(value)

    def setDiffusionCoefficient(self, value):
        self._diffusionCoefficient = value

    def getDepField(self):
        return self.depField

    def correctBCs(self):

        for loc, typeValueTupel in self.depField.boundary.items():
            (type, value) = typeValueTupel
            if type == 'zeroGradient':
                self.setVonNeumann(loc)
            else:
                self.setDerichlet(loc, value)

    # # def updateLinSystem(self):
    # #     self._linSystem.reset(shape=self._depFieldShape)
    # #     self._centreMatrixCoeffs = self._linSystem.update(F=self._convFluxes, D=self._diffFluxes, Sc=self._sourceField_c, Sp=self._sourceField_p)
    #
    # def getCentreMatrixCoeffs(self):
    #     return self._centreMatrixCoeffs

    def solve(self):
        return self._linSystem.solve()

    def setVonNeumann(self, loaction):
        pass

    def setDerichlet(self, location, value):
        pass
