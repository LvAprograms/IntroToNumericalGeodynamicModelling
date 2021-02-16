from numpy import zeros, ones, linspace
from scipy import linalg
from matplotlib import pyplot as plt


# Model parameters
Lx = 1e6 # model length in m
nx = 101 # amount of nodal points
x = linspace(0, Lx, nx)
dx = Lx / (nx -1) # grid step
"""
Phi0 ---dPhiA ---ddPhi1 ---
x0 ----- x1 ----- x2 ----- ... ----- x_nx-1
     xA ----- xB ----- xC ----- ...
"""

phi = zeros((1, nx))    # potential
dphi = zeros((1,nx-2))    # first derivative of potential
ddphi = ones((1, nx-4))  # second derivative of potential
# boundary conditions:
phi[0] = 0
phi[-1] = 0

L = zeros((nx, nx))
r = ones((nx, 1))
r[0] = 0
r[-1] = 0
L[0,0] = 1
L[-1, -1] = 1

# the coefficients are what's in front of Phi1, Phi2 etc. Since dx = x[i] - x[i-1], it comes down to this:
# 2 / (x2 - x1) / (x3 - x1) = 2 / dx / (2 * dx) = 2 / (dx * 2 * dx) = 2 / (2 * dx ** 2) = 1 / dx ** 2
for i in range(1, nx-1):
    # fill left and righthand side?
    L[i, i-1] = 1 / dx ** 2
    L[i, i] =  -2 / dx ** 2
    L[i, i+1] = 1 / dx ** 2
    r[i] = 1


# solve directly
S = linalg.solve(L, r)

# plot results
plt.figure()
plt.plot(x, S)
plt.xlim([0, Lx])
plt.ylim([min(S), 0])
plt.xlabel("Distance [m]")
plt.ylabel("Gravitational potential [?]")
plt.grid(b=True)
plt.show()