import Fields

def scalarDiffusion(mesh, field, diffCoeff):
    ### returning coefficient vectors for subsequent coefficient matrix assembly and coeff vector

    faceFluxes = Fields.parameterFaceField(mesh=mesh)
    gamma = Fields.parameterFaceField(mesh=mesh, value=diffCoeff)

    # a mesh method
    faceArea = Fields.parameterFaceField(mesh=mesh, value=1.0)

    # constant boundary values, should be read from field parameter
    Tw = 100
    Te = 500
    Tn = 2
    Ts = 1

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

    s_u = Fields.parameterCellField(mesh=mesh, value=0)

    a_p = -(a_e + a_w + a_n + a_s)

    # defining neumann BC at north:
    a_p.bn += a_n.bn
    s_u.bn += 0

    a_n.bn = 0

    #defining vNeumann BC at south:
    a_p.bs += a_s.bs
    s_u.bs += 0

    a_s.bs = 0

    # defining Dirichlet BC at east/west boundaries:
    # should be read from field parameter

    s_u.be -= a_e.be * Te
    a_e.be = 0

    s_u.bw -= a_w.bw * Tw
    a_w.bw = 0


    return a_e.raw, a_w.raw, a_n.raw, a_s.raw, a_p.raw, s_u.raw
