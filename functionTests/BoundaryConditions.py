

class staggeredBC:

    @staticmethod
    def derichlet(loc, value, **kwargs):
        D = kwargs.get('D')
        F = kwargs.get('F')
        Sc = kwargs.get('Sc')
        Sp = kwargs.get('Sp')

        if loc == 'left':
            Sc.bw += 1e15 * value
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

    @staticmethod
    def vonNeumann(loc, **kwargs):
        D = kwargs.get('D')
        F = kwargs.get('F')
        Sc = kwargs.get('Sc')
        Sp = kwargs.get('Sp')

        # I just set these fluxes to zero. But I should set them to their intenal neighbors value
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


class scalarBC:

    @staticmethod
    def derichlet(loc, value, **kwargs):
        D = kwargs.get('D')
        F = kwargs.get('F')
        Sc = kwargs.get('Sc')
        Sp = kwargs.get('Sp')

        if loc == 'left':
            Sc.bw += (2.0 * D.u.bw + F.u.bw) * value
            Sp.bw += 2.0 * D.u.bw + F.u.bw
            D.u.bw = 0.0
            F.u.bw = 0.0
        elif loc == 'right':
            Sc.be += (2.0 * D.u.be - F.u.be) * value
            Sp.be += 2.0 * D.u.be - F.u.be
            D.u.be = 0.0
            F.u.be = 0.0
        elif loc == 'top':
            Sc.bn += (2.0 * D.v.bn - F.v.bn) * value
            Sp.bn += 2.0 * D.v.bn - F.v.bn
            D.v.bn = 0.0
            F.v.bn = 0.0
        elif loc == 'bottom':
            Sc.bs += (2.0 * D.v.bs + F.v.bs) * value
            Sp.bs += 2.0 * D.u.bs + F.v.bs
            D.v.bs = 0.0
            F.v.bs = 0.0

    # why are convective fluxes not affected?
    @staticmethod
    def vonNeumann(loc, **kwargs):
        D = kwargs.get('D')
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
