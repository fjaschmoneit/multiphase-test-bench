'''
This is a validation test for the correct implementation of an implicit time stepping scheme.
The problem follows the example in Versteeg/Malalasekera, ch. 8, example 8.2.

A 2cm long metal plate is held 200 deg C.
The temperature on its right edge drops suddenly to 0 deg C.
The left side is insulated, a zero gradient BC is assumed here.

Because the temperature on the left side in not held constant, this problem is transient.
The equated temperature profiles compare well with the analytical solutions.
'''

import Manager as Odin
import TransportModels
import matplotlib.pyplot as plt
import numpy as np

# geometric parameters
LenX = 0.02
nbCells = 5
resolution = nbCells/LenX

# thermal conductivity [ W / (m K) ]
k = 10

# thermal capacity times density [ J / (m3 K) ]
rhoc = 1e7

# initial temperature
T0 = 200

# time step [s]
dt = 2
endTime = 160

# set cross sectional area to 1
geom = Odin.createGeometry( 'rectangle', [LenX, LenX/5] )
mesh = Odin.createMesh( geom, res=resolution )

flowModels = {
    'T' : TransportModels.scalarTransport
}

Odin.initialize(flowModels, mesh, geom, passiveFields={'u' : 'faces_u','v' : 'faces_v'} )

T = Odin.getField('T')

T.govModel.setDiffusionCoefficient(k)
T.govModel.setHeatCapacity(rhoc)

T.data.fill(T0)
Odin.defineBoundaryCondition(field=T, boundaryName='top', type='zeroGradient' )
Odin.defineBoundaryCondition(field=T, boundaryName='bottom', type='zeroGradient' )
Odin.defineBoundaryCondition(field=T, boundaryName='left', type='zeroGradient' )
Odin.defineBoundaryCondition(field=T, boundaryName='right', type='fixedValue', value=0 )

# initiating time loop:
currentTime = 0
while currentTime <= endTime:
    T.data = Odin.transientSolve(T, dt, method='implicit')
    if currentTime%20 == 0:
        plt.plot(np.linspace(0.5*LenX/nbCells,LenX*(1-0.5/nbCells),nbCells), T.data[0], '-o', label=f'{currentTime} s')
    currentTime += dt

# post-processing:
print(T.data)
ax = plt.gca()
ax.set(xlim=(0, LenX), ylim=(0, 1.1*T0))
ax.set_xticks(np.linspace(0,LenX,nbCells+1))
plt.legend()
plt.grid()
plt.show()
