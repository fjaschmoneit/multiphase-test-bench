import sys
sys.path.append('../src/MultiphaseTestBench')
import Manager

LenX = 5
LenY = 1
resolution = 1

T_l = 0.1
T_r = 0.5


mesh = Manager.createMesh( Manager.createGeometry( 'rectangle', [LenX, LenY] ), res=resolution )

T = Manager.newField(type='scalar',
                  transportModels={Manager.TransportModels_dev.Diffusion})



# assembling simulation instance:

#
# Manager.newSim
#
# Manager.objReg.MESH =
# Manager.objReg.FIELDS =



print("moin")