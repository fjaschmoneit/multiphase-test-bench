import LinearEquationSystems
import Mesh

mesh = Mesh.cartesian2D(1,0.1,10)
Sys = LinearEquationSystems.linearSystem(mesh)
Sys.update()

x = Sys.solve()
print(x)