from fieldAccess import*

class scalarBC:

    @staticmethod
    def derichlet(loc, value, transportInstance):

        if loc == 'left':
            westFlux = transportInstance.a_w[boundary_west]
            transportInstance.a_p[boundary_west] += 2*westFlux
            transportInstance.sourceField_c[boundary_west] += 2*westFlux*value
            transportInstance.a_w[boundary_west] = 0.0

        elif loc == 'right':
            eastFlux = transportInstance.a_e[boundary_east]
            transportInstance.a_p[boundary_east] += 2*eastFlux
            transportInstance.sourceField_c[boundary_east] += 2*eastFlux*value
            transportInstance.a_e[boundary_east] = 0.0

        elif loc == 'top':
            northFlux = transportInstance.a_n[boundary_north]
            transportInstance.a_p[boundary_north] += 2*northFlux
            transportInstance.sourceField_c[boundary_north] += 2*northFlux*value
            transportInstance.a_n[boundary_north] = 0.0

        elif loc == 'bottom':
            southFlux = transportInstance.a_s[boundary_south]
            transportInstance.a_p[boundary_south] += 2*southFlux
            transportInstance.sourceField_c[boundary_south] += 2*southFlux*value
            transportInstance.a_s[boundary_south] = 0.0

    # why are convective fluxes not affected?
    @staticmethod
    def vonNeumann(loc, transportInstance):
        if loc == 'left':
            transportInstance.a_w[boundary_west].fill(0)
        elif loc == 'right':
            transportInstance.a_e[boundary_east].fill(0)
        elif loc == 'top':
            transportInstance.a_n[boundary_north].fill(0)
        elif loc == 'bottom':
            transportInstance.a_s[boundary_south].fill(0)
