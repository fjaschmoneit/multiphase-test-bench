from fieldAccess import*

class scalarBC:

    @staticmethod
    def derichlet(loc, value, transportInstance):

        if loc == 'left':
            westFlux = transportInstance.a_w[boundary_west]
            transportInstance.a_p[boundary_west] += 2*westFlux
            transportInstance._sourceField_c[boundary_west] += 2*westFlux*value
            transportInstance.a_w[boundary_west] = 0.0

        elif loc == 'right':
            eastFlux = transportInstance.a_e[boundary_east]
            transportInstance.a_p[boundary_east] += 2*eastFlux
            transportInstance._sourceField_c[boundary_east] += 2*eastFlux*value
            transportInstance.a_e[boundary_east] = 0.0

        elif loc == 'top':
            northFlux = transportInstance.a_n[boundary_north]
            transportInstance.a_p[boundary_north] += 2*northFlux
            transportInstance._sourceField_c[boundary_north] += 2*northFlux*value
            transportInstance.a_n[boundary_north] = 0.0

        elif loc == 'bottom':
            southFlux = transportInstance.a_s[boundary_south]
            transportInstance.a_p[boundary_south] += 2*southFlux
            transportInstance._sourceField_c[boundary_south] += 2*southFlux*value
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


class staggeredBC:

    # for u-cell only
    @staticmethod
    #def derichlet(loc, value, transportInstance):
        #
        # # fixing normal velocity:
        # if loc == 'left':
        #     transportInstance.a_p[boundary_west] += 1e15
        #     transportInstance._sourceField_c[boundary_west] += 1e15 * value
        #     transportInstance.a_w[boundary_west] = 0.0
        # elif loc == 'right':
        #     transportInstance.a_p[boundary_east] += 1e15
        #     transportInstance._sourceField_c[boundary_east] += 1e15 * value
        #     transportInstance.a_e[boundary_east] = 0.0



        # elif loc == 'top':
        #     northFlux = transportInstance.a_n[boundary_north]
        #
        #     #shearSource =
        #     transportInstance.a_p[boundary_north] += 2*northFlux
        #     transportInstance._sourceField_c[boundary_north] += 2*northFlux*value
        #     transportInstance.a_n[boundary_north] = 0.0
        #
        # elif loc == 'bottom':
        #     southFlux = transportInstance.a_s[boundary_south]
        #     transportInstance.a_p[boundary_south] += 2*southFlux
        #     transportInstance._sourceField_c[boundary_south] += 2*southFlux*value
        #     transportInstance.a_s[boundary_south] = 0.0

        #
        # elif loc == 'top':
        #     transportInstance.a_p[boundary_north] += 1e15
        #     transportInstance._sourceField_c[boundary_north] += 1e15 * value
        #     transportInstance.a_n[boundary_north] = 0.0
        # elif loc == 'bottom':
        #     transportInstance.a_p[boundary_south] += 1e15
        #     transportInstance._sourceField_c[boundary_south] += 1e15 * value
        #     transportInstance.a_s[boundary_south] = 0.0

    @staticmethod
    def vonNeumann(loc, transportInstance):
        return
        if loc == 'left':
            transportInstance.a_w[boundary_west] = -transportInstance.a_w[boundary_nb1_west]
        elif loc == 'right':
#            transportInstance.a_w[boundary_west] = transportInstance.a_w[boundary_nb1_west]
            transportInstance.a_e[boundary_east] = 0.5*transportInstance.a_w[boundary_east]

        # elif loc == 'top':
        #     transportInstance.a_n[boundary_north] = -transportInstance.a_n[boundary_nb1_north]
        # elif loc == 'bottom':
        #     transportInstance.a_s[boundary_south] = -transportInstance.a_s[boundary_nb1_south]

        # D = kwargs.get('D')
        # F = kwargs.get('F')
        # Sc = kwargs.get('Sc')
        # Sp = kwargs.get('Sp')
        #
        # # I just set these fluxes to zero. But I should set them to their intenal neighbors value
        # if loc == 'left':
        #     D.u.bw = 0.0
        #     F.u.bw = 0.0
        # elif loc == 'right':
        #     D.u.be = 0.0
        #     F.u.be = 0.0
        # elif loc == 'top':
        #     D.v.bn = 0.0
        #     F.v.bn = 0.0
        # elif loc == 'bottom':
        #     D.v.bs = 0.0
        #     F.v.bs = 0.0
        #

