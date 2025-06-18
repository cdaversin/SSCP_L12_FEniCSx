import ufl
from ufl import TestFunction, TrialFunction, dot, inner, grad, dx, sqrt

from dolfinx.mesh import create_unit_square, CellType
from dolfinx.fem import Constant, functionspace, Function, locate_dofs_geometrical, dirichletbc, Expression
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile

from mpi4py import MPI
import matplotlib.pyplot as plt
import numpy as np
import pyvista

# Unit square mesh
mesh = create_unit_square(comm=MPI.COMM_WORLD, nx=20, ny=20, cell_type=CellType.triangle)
# Lagrange (CG1) function space
W = functionspace(mesh, ("Lagrange", 1))

# kappa / mu
kappa_over_mu = Constant(mesh, 1.0)
# porosity, ranging from 0 to 1
phi = Constant(mesh, 0.1)
# Source term
S = Constant(mesh, 0.0)

# Variational formulation
v = TestFunction(W)
p = TrialFunction(W)

a = kappa_over_mu * dot(grad(p), grad(v)) * dx
L = S *  v * dx

# Dirichlet boundary conditions
value_left = Constant(mesh, 1.0)
value_right = Constant(mesh, 0.0)

# Imposing Dirichlet BC to the left boundary node
dofs_left = locate_dofs_geometrical(W, lambda x: np.isclose(x[0], 0))
bc_left = dirichletbc(value_left, dofs_left, W)
# Imposing Dirichlet BC to the right boundary node
dofs_right = locate_dofs_geometrical(W, lambda x: np.isclose(x[0], 1))
bc_right = dirichletbc(value_right, dofs_right, W)

bcs = [bc_left, bc_right]

# Solve
problem = LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
ph = problem.solve()

# Plot the velocity (post-processing)
# Velocity is a vector - need a vectorial function space
gdim = mesh.geometry.dim
W_v = functionspace(mesh, ("DG", 0, (gdim,)))
W_s = functionspace(mesh, ("DG", 0))
vh = Function(W_v)
v_expr = Expression(-kappa_over_mu * grad(ph) / phi, W_v.element.interpolation_points())
vh.interpolate(v_expr)

# Magnitude
magnitude_expr = Expression(sqrt(inner(vh, vh)), W_s.element.interpolation_points())
vh_magnitude = Function(W_s)
vh_magnitude.interpolate(magnitude_expr)

# Plot - pyvista
from dolfinx import plot

# Create the grid to plot the mesh
mesh.topology.create_connectivity(2, 2)
topology, cell_types, geometry = plot.vtk_mesh(mesh, 2)
p_grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
# Plot solution on the mesh
p_grid.point_data["p"] = ph.x.array
p_grid.set_active_scalars("p")
p_plotter = pyvista.Plotter()
p_plotter.add_mesh(p_grid, show_edges=True)
p_plotter.view_xy()
p_plotter.save_graphic('Darcy_p_2D.pdf')

v_grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
# Plot solution on the mesh
v_grid.cell_data["v"] = vh_magnitude.x.array
v_grid.set_active_scalars("v")
v_plotter = pyvista.Plotter()
v_plotter.add_mesh(v_grid, show_edges=True)
v_plotter.view_xy()
v_plotter.save_graphic('Darcy_v_magnitude_2D.pdf')


# Plot velocity (xdmf) - visualize with paraview
with XDMFFile(MPI.COMM_WORLD, "Darcy_velocity.xdmf", "w") as xdmf_file:
    # Write the mesh (only needs to be done once per file)
    xdmf_file.write_mesh(mesh)
    # Write scalar function
    xdmf_file.write_function(vh)
