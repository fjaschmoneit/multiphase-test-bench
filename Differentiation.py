from fieldAccess import *


def grad_u(field, mesh):
    rCellDist_u = mesh.calcInvCellDistance('west')
    return (field[east] - field[west])*rCellDist_u[internal_u]


def grad_v(field, mesh):
    rCellDist_v = mesh.calcInvCellDistance('north')
    return (field[north] - field[south])*rCellDist_v[internal_v]



#
# def grad(cellField, mesh):
#     # differentialtion of cellField
#     # returns a parameterFaceField
#
#     ff = Fields.parameterFaceField(mesh=mesh, value=0)
#     invCellDist = mesh.invCellDist
#
#     d_x = (cellField.raw[:,1:] - cellField.raw[:,:-1])*invCellDist.internalEntries_EW
#     d_y = (cellField.raw[1:,:] - cellField.raw[:-1,:])*invCellDist.internalEntries_NS
#
#     ff.internalEntries_EW = d_x
#     ff.internalEntries_NS = d_y
#
#     return ff
#
#
