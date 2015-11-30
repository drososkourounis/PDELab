"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Simplest example of computation and visualization with FEniCS.

-Laplace(u) = f on the unit square.
u = u0 on the boundary.
u0 = u = 1 + x^2 + 2y^2, f = -6.
"""

from dolfin import *

# Create mesh and define function space
#mesh = UnitSquare(6, 4)
info(parameters, True)
set_log_level(DEBUG)
list_krylov_solver_preconditioners()
mesh = UnitCubeMesh(50, 50, 50)
V = FunctionSpace(mesh, 'Lagrange', 1)

# Define boundary conditions
u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]+x[2]*x[2]')

def u0_boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u0, u0_boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0) # edo prepei na einai delta function
a = inner(grad(u), grad(v))*dx

L = f*v*dx

# Compute solution
u = Function(V)
problem = LinearVariationalProblem(a, L, u, bc)
solver  = LinearVariationalSolver(problem)
solver.parameters['linear_solver']  = "cg"
solver.parameters['preconditioner'] = "petsc_amg"
gmres_params = solver.parameters["krylov_solver"]
gmres_params['relative_tolerance'] = 1.0e-6
gmres_params['maximum_iterations'] = 1000
gmres_params.gmres['restart'] = 1000
gmres_params['monitor_convergence'] = True

precond_params = gmres_params['preconditioner']
precond_params.schwarz['overlap'] = 1

solver.solve()

# Plot solution and mesh
#plot(u)
#plot(mesh)

# Dump solution to file in VTK format
#file = File('poisson.pvd')
#file << u

# Hold plot
#interactive()
