# Here we put some useful classes and functions for the PDE-regularized prior projects
# Chao Zhang, 31.01.2024

import dolfin as df
import numpy as np
from numpy import linalg as LA

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
        return self.dudf.T@y
