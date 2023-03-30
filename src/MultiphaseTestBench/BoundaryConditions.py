from fieldAccess import*
import Interpolation


LARGE = 1e20


class scalarBC:

    # @staticmethod
    # def derichlet_(transportInstance, **argDict):
    #     # second order accurate
    #     direction = argDict.get('direction')
    #     value = argDict.get('value')
    #
    #     (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
    #     a_dir = transportInstance.getCoefficientsInDirection(direction)
    #     a_opp_dir = transportInstance.getCoefficientsInOppositeDirection(direction)
    #
    #     D = transportInstance.diffusionCoefficient
    #     fA = transportInstance.mesh.calcFaceArea(direction=direction)[boundary_dir]
    #     invCellDist = transportInstance.mesh.calcInvCellDistance(direction=direction)[boundary_dir]
    #     b = D*fA*invCellDist/3.0
    #
    #     a_dir[boundary_dir] = 0
    #     transportInstance.a_p[boundary_dir] += 8.0*b
    #     a_opp_dir[boundary_dir] += b
    #     transportInstance.sourceField_c[boundary_dir] += 8.0*b*value + transportInstance.u.data[boundary_dir]

    # @staticmethod
    # def vonNeumann(transportInstance, **argDict):
    #     direction = argDict.get('direction')
    #
    #     (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
    #     a_dir = transportInstance.getCoefficientsInDirection(direction)
    #     a_opp_dir = transportInstance.getCoefficientsInOppositeDirection(direction)
    #     fA = transportInstance.mesh.calcFaceArea(direction=direction)[boundary_dir]
    #
    #     a_dir[boundary_dir] = 0.0
    #     a_opp_dir[boundary_dir] += 0.5*transportInstance.v.data[boundary_dir]*fA[boundary_dir]

        #
        #
        #
        #
        # a_p = transportInstance.getCentreMatrixCoeffs()
        # S_c = transportInstance.sourceField_c
        #
        # ghostflux = a_dir[boundary_dir]
        # a_p[boundary_dir] += 2.0*ghostflux
        # S_c[boundary_dir] += 2.0*ghostflux*value
        # a_dir[boundary_dir] = 0.0


    @staticmethod
    def derichlet(transportInstance, **argDict):
        # first order accurate
        direction = argDict.get('direction')
        value = argDict.get('value')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_p = transportInstance.getCentreMatrixCoeffs()
        S_c = transportInstance.sourceField_c

        ghostflux = a_dir[boundary_dir]
        a_p[boundary_dir] += 2.0*ghostflux
        S_c[boundary_dir] += 2.0*ghostflux*value
        a_dir[boundary_dir] = 0.0

    @staticmethod
    def vonNeumann(transportInstance, **argDict):
        direction = argDict.get('direction')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_dir[boundary_dir] = 0.0

    @staticmethod
    def symmetry(transportInstance, **argDict):
        direction = argDict.get('direction')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_dir[boundary_dir] = 0.0
        #transportInstance.a_p[boundary_dir] +=


class staggered_u:

    @staticmethod
    def vonNeumann(transportInstance, **argDict):
        direction = argDict.get('direction')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        # a_dir[boundary_dir] = 0.0

#        transportInstance.a_p[boundary_dir] += a_dir[boundary_nb1_dir]

        #not so good:
        # a_opposite = self.getCoefficientsInDirection(opposite(direction))
        # a_opposite[boundary_dir] = 1e20
        a_dir[boundary_dir] = 0.0
        #a_dir[boundary_dir] = self.a_p[boundary_dir]

        #self.a_w[boundary_dir] = 1e20

    @staticmethod
    def derichlet(transportInstance, **argDict):
        direction = argDict.get('direction')
        value = argDict.get('value')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_p = transportInstance.getCentreMatrixCoeffs()
        S_c = transportInstance.sourceField_c

        if direction == 'west' or direction == 'east':
            a_p[boundary_dir] += LARGE
            S_c[boundary_dir] += LARGE * value
            a_dir[boundary_dir] = 0.0

        elif direction == 'north' or direction == 'south':

            faceAreas_u, faceAreas_v = transportInstance.mesh.calcFaceAreas()
            rCellDist_u, rCellDist_v = transportInstance.mesh.calcInvCellDistances()

            invCellWallDist_A = 2.0 * rCellDist_v[boundary_dir]*faceAreas_v[boundary_dir]
            nodeFluxes = Interpolation.toStaggered(invCellWallDist_A, 'u')

            a_p[boundary_dir] += transportInstance.diffusionCoefficient*nodeFluxes
            S_c[boundary_dir] += transportInstance.diffusionCoefficient*nodeFluxes*value
            a_dir[boundary_dir] = 0.0

class staggered_v:

    @staticmethod
    def vonNeumann(transportInstance, **argDict):
        direction = argDict.get('direction')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_dir[boundary_dir] = 0.0

    @staticmethod
    def derichlet(transportInstance, **argDict):
        direction = argDict.get('direction')
        value = argDict.get('value')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_p = transportInstance.getCentreMatrixCoeffs()
        S_c = transportInstance.sourceField_c

        if direction == 'north' or direction == 'south':
            a_p[boundary_dir] += LARGE
            S_c[boundary_dir] += LARGE * value
            a_dir[boundary_dir] = 0.0

        elif direction == 'east' or direction == 'west':

            faceAreas_u, faceAreas_v = transportInstance.mesh.calcFaceAreas()
            rCellDist_u, rCellDist_v = transportInstance.mesh.calcInvCellDistances()

            invCellWallDist_A = 2.0 * rCellDist_u[boundary_dir] * faceAreas_u[boundary_dir]
            nodeFluxes = Interpolation.toStaggered(invCellWallDist_A, 'v')

            a_p[boundary_dir] += transportInstance.diffusionCoefficient * nodeFluxes
            S_c[boundary_dir] += transportInstance.diffusionCoefficient * nodeFluxes * value
            a_dir[boundary_dir] = 0.0

class pressure:

    @staticmethod
    def freeFlow(transportInstance, **argDict):
        direction = argDict.get('direction')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_dir[boundary_dir] = 0.0

    @staticmethod
    def totalPressure(transportInstance, **argDict):
        pass

    @staticmethod
    def derichlet(transportInstance, **argDict):
        direction = argDict.get('direction')
        value = argDict.get('value')

        (boundary_dir, boundary_nb1_dir) = fieldSlice(direction)
        a_dir = transportInstance.getCoefficientsInDirection(direction)
        a_p = transportInstance.getCentreMatrixCoeffs()
        S_c = transportInstance.sourceField_c

        a_p[boundary_dir] += LARGE
        S_c[boundary_dir] += LARGE * value
        a_dir[boundary_dir] = 0.0

