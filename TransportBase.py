import Fields

class transportBase:

    def __init__(self, mesh, fieldCreator, fieldReg, linSystem, depFieldShape):
        self.mesh = mesh
        self.linSystem = linSystem
        self.fieldReg = fieldReg
        self.fc = fieldCreator
        self.depFieldShape = depFieldShape

        self.depField = self.fc.newField(shape=self.depFieldShape, value=0.0)  # why is this not a simple data field?
        self.sourceField_c = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.a_e = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.a_w = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.a_s = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.a_n = Fields.newDataField(shape=self.depFieldShape, value=0.0)
        self.a_p = Fields.newDataField(shape=self.depFieldShape, value=0.0)

    def reset(self):
#        self.depField.data.fill(0.0)
        self.sourceField_c.fill(0.0)
        self.a_e.fill(0.0)
        self.a_w.fill(0.0)
        self.a_s.fill(0.0)
        self.a_n.fill(0.0)
        self.a_p.fill(0.0)

    def setSourceField(self, value):
        self.sourceField_c.fill(value)

    def setDiffusionCoefficient(self, value):
        self.diffusionCoefficient = value

    def getCentreMatrixCoeffs(self):
        return self.a_p

    def getDepField(self):
        return self.depField

    def correctBCs(self):

        for loc, typeValueTupel in self.depField.boundary.items():
            (type, value) = typeValueTupel
            if type == 'zeroGradient':
                self.setVonNeumann(loc)
            else:
                self.setDerichlet(loc, value)

    def updateLinSystem(self):
        self.linSystem.reset(shape=self.depFieldShape)
        self.linSystem.set_e_coeffs(self.a_e)
        self.linSystem.set_w_coeffs(self.a_w)
        self.linSystem.set_n_coeffs(self.a_n)
        self.linSystem.set_s_coeffs(self.a_s)
        self.linSystem.set_p_coeffs(self.a_p)
        self.linSystem.set_b(self.sourceField_c)

    # def getCentreMatrixCoeffs(self):
    #     return self.centreMatrixCoeffs

    def solve(self):
        return self.linSystem.solve()

    def setVonNeumann(self, loaction):
        pass

    def setDerichlet(self, location, value):
        pass
