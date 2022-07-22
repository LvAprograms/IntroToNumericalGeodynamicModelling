from numpy import zeros, ones, linspace, meshgrid
from scipy import linalg
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

# Model parameters
Lx = 1000 # model size in horizontal direction
Ly = 1500 # Model size in vertical direction
nx = 31 # grid points in x
ny = 41 # grid points in y
dx = Lx / (nx - 1)
dy = Ly / (ny - 1)
x = linspace(0, Lx, nx)
y = linspace(0, Ly, ny)
N = nx * ny # dimension of global matrix L

L = zeros((N, N))
R = zeros((N, 1)) # right hand side vector

"""
j-1             j           j+1
 ------------Phi[i-1,j]----|-------- i-1
 |              |          |
 |              |          |
 Phi[i,j-1]--Phi[i,j]---Phi[i,j+1] -- i
 |              |          | 
 |              |          |
 ------------Phi[i+1,j]--------------i+1

"""
for i in range(ny):
    for j in range(nx):
        k = (i) * nx + j # other way around because python arrays are row major while Matlab arrays are column major
        if i == 0 or i == ny-1 or j == 0 or j == nx-1:
            # you're at the boundary
            L[k,k] = 1 
            R[k] = 0
        else:
            L[k, k - 1] = 1 / dy ** 2                   # get coefficient for Phi[i-1, j]
            L[k, k - nx]  = 1 / dx ** 2                 # get coefficient for Phi[i, j-1]
            L[k, k]      = -2 / dx ** 2 - 2 / dy ** 2   # get coefficient for Phi[i,j]
            L[k, k + nx] = 1 / dx ** 2                  # get coefficient for Phi[i, j + 1]
            L[k, k + 1]  = 1 / dy ** 2                  # get coefficient for Phi[i+1, j]
            R[k]         = 1

# solve directly
S = linalg.solve(L, R)

xx, yy = meshgrid(x, y)
fig = plt.figure(figsize=(30,10))
ax0 = fig.add_subplot(1,3,1, projection='3d')
ax1 = fig.add_subplot(1,3,2, projection='3d')
ax2 = fig.add_subplot(1,3,3, projection='3d')

cf = ax0.plot_surface(xx, yy, S.reshape((ny, nx)), cmap='hot')
fig.colorbar(cf, ax=ax0, orientation='horizontal', label="Gravitational potential")

for ax in [ax0, ax1, ax2]:
    ax.set_xlim([0,1000])
    ax.set_ylim([0,1500])
    ax.set_zlim([-100000, 0])
# plt.show()
ax0.set_title("Direct")
# Exercise 3.3 Gauss-Seidel iterations
Phi      = zeros((ny, nx)) # Estimate unknowns to be zero
Snew     = zeros((ny, nx))     
deltaR   = zeros((ny, nx))
thetaGS  = 1.5
thetaJK  = 1.0
plotinterval = 10
niter = 300
R = ones((ny, nx))
R[0,:], R[:,0], R[-1,:], R[:,-1] = 0, 0, 0, 0
for iter in range(niter):
    for i in range(1, ny-1):
        for j in range(1, nx-1):
            deltaR[i,j] = R[i,j] - ((Phi[i, j-1] - 2 * Phi[i,j] + Phi[i, j+1]) / dx ** 2 + (Phi[i-1,j] - 2 * Phi[i, j]+ Phi[i+1,j]) / dy ** 2)
            # For Gauss-Seidel iterations, the new values are immediately assigned to the solution estimate
            Phi[i,j] += thetaGS * deltaR[i,j] / (-2 / dx**2 - 2 / dy ** 2)
    if (iter+1) % plotinterval == 0:
        print("plotting results after iteration {}".format(iter+1))
        ax1.plot_surface(xx,yy,Phi, cmap='hot')
        resid = ax1.plot_surface(xx,yy,deltaR, cmap='cool')
fig.colorbar(resid, orientation='horizontal', ax=ax1, label=" Residual")
ax1.set_title("300 Gauss-Seidel iterations ")


# Exercise 3.4 Jacobi iterations
Phi      = zeros((ny, nx)) # Estimate unknowns to be zero
Snew     = zeros((ny, nx))     
deltaR   = zeros((ny, nx))
thetaJK  = 1.0
plotinterval = 10
niter = 300
R = ones((ny, nx))
R[0,:], R[:,0], R[-1,:], R[:,-1] = 0, 0, 0, 0
for iter in range(niter):
    for i in range(1, ny-1):
        for j in range(1, nx-1):
            deltaR[i,j] = R[i,j] - ((Phi[i, j-1] - 2 * Phi[i,j] + Phi[i, j+1]) / dx ** 2 + (Phi[i-1,j] - 2 * Phi[i, j]+ Phi[i+1,j]) / dy ** 2)
    # Jacobi iterations: the new values are assigned after the entire solution is estimated
    Phi += thetaJK * deltaR / (-2 / dx**2 - 2 / dy ** 2)
    
    if (iter+1) % plotinterval == 0:
        print("plotting results after iteration {}".format(iter+1))
        ax2.plot_surface(xx,yy,Phi, cmap='hot')
        resid = ax2.plot_surface(xx,yy,deltaR, cmap='cool')
        # plt.pause(0.001)
fig.colorbar(resid, orientation='horizontal', ax=ax2, label=" Residual")
ax2.set_title("300 Jacobi iterations")
plt.show()