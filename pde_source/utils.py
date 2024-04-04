# Here we put some useful classes and functions for the PDE-regularized prior projects
# Chao Zhang, 31.01.2024

import dolfin as df
import numpy as np

class Poisson:
    def __init__(self, nx):
        self.nx = nx
        self.mesh = df.UnitIntervalMesh(self.nx)
        self.V = df.FunctionSpace(self.mesh, "Lagrange", 1)
        def boundary(x):
            return x[0] < df.DOLFIN_EPS or x[0] > 1.0 - df.DOLFIN_EPS
        u0 = df.Constant(0.0)
        self.bc = df.DirichletBC(self.V, u0, boundary)

        self.u = df.TrialFunction(self.V)
        self.v = df.TestFunction(self.V)
        self.f_exact = df.Expression("1.0*(x[0]<0.6)*(x[0]>0.4)+0.0", degree=1)
        self.f_exact = df.interpolate(self.f_exact, self.V)

        self.a = df.inner(df.grad(self.u), df.grad(self.v))*df.dx

        self.A = df.assemble(self.a)

        self.bc.apply(self.A)

        self.f = df.Function(self.V)

        self.sol = df.Function(self.V)

        # related to gradient
        self.M = df.assemble(df.inner(self.u, self.v)*df.dx)
        self.bc.apply(self.M)

        self.AA = self.A.array()
        self.MM = self.M.array()

        self.MM[0,:] = 0
        self.MM[-1,:] = 0

        self.AA_inv = np.linalg.inv(self.AA)

        self.dudf = self.AA_inv@self.MM

    def solve_exact(self):
        L = self.f_exact*self.v*df.dx

        b = df.assemble(L)

        self.bc.apply(b)

        df.solve(self.A, self.sol.vector(), b)

        return self.sol.vector()[:]

    def forward(self, f_numpy):
        self.f.vector()[:] = f_numpy
        L = self.f*self.v*df.dx

        b = df.assemble(L)

        self.bc.apply(b)

        df.solve(self.A, self.sol.vector(), b)

        return self.sol.vector()[:].T

    def adjoint(self, y):
        return self.dudf@y

class Poisson2D:
    def __init__(self, nx, ny):
        self.nx = nx
        self.ny = ny
        self.mesh = df.UnitSquareMesh(self.nx, self.ny)
        self.V = df.FunctionSpace(self.mesh, "Lagrange", 1)
        def left(x, on_boundary):
            return on_boundary and x[0] < df.DOLFIN_EPS
        def right(x, on_boundary):
            return on_boundary and x[0] > 1.0 - df.DOLFIN_EPS
        def bottom(x, on_boundary):
            return on_boundary and x[1] < df.DOLFIN_EPS
        def top(x, on_boundary):
            return on_boundary and x[1] > 1.0 - df.DOLFIN_EPS
        def boundary(x):
            return x[0] < df.DOLFIN_EPS or x[0] > 1.0 - df.DOLFIN_EPS
        ub_left = df.Expression("1.0*(x[1]<0.6)*(x[1]>0.4)+0.0", degree=1)
        ub_right = df.Constant(0.0)
        ub_bottom = df.Constant(0.0)
        ub_top = df.Constant(0.0)
        self.bc_left = df.DirichletBC(self.V, ub_left, left)
        self.bc_right = df.DirichletBC(self.V, ub_right, right)
        self.bc_bottom = df.DirichletBC(self.V, ub_bottom, bottom)
        self.bc_top = df.DirichletBC(self.V, ub_top, top)
        self.bcs = [self.bc_left, self.bc_right, self.bc_bottom, self.bc_top]

        self.left_bc_dofs = np.array(list(self.bc_left.get_boundary_values().keys()))

        self.u = df.TrialFunction(self.V)
        self.v = df.TestFunction(self.V)
        self.f = df.Expression("0.0", degree=1)
        self.f = df.interpolate(self.f, self.V)

        self.a = df.inner(df.grad(self.u), df.grad(self.v))*df.dx

        self.A = df.assemble(self.a)

        for bc in self.bcs:
            bc.apply(self.A)

        # self.f = df.Function(self.V)

        self.sol = df.Function(self.V)

        self.AA = self.A.array()

        self.AA_inv = np.linalg.inv(self.AA)

    def solve_exact(self):
        L = self.f*self.v*df.dx

        b = df.assemble(L)

        for bc in self.bcs:
            bc.apply(b)

        df.solve(self.A, self.sol.vector(), b)

        return self.sol.vector()[:]

    def forward(self, bc_numpy):
        L = self.f*self.v*df.dx

        b = df.assemble(L)

        for bc in self.bcs:
            bc.apply(b)
        
        b[self.left_bc_dofs] = bc_numpy

        df.solve(self.A, self.sol.vector(), b)

        return self.sol.vector()[:].T

    def adjoint(self, y):
        return self.AA_inv[:,self.left_bc_dofs].T@y
        # return self.dudf@y