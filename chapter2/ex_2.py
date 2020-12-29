from numpy import zeros, arange, meshgrid, mean, min, max
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from chapter2.gibbsfree import Gibbs_m


## Exercise 2.2 (2.1 is analytical)
T0 = 100    # deg K
T1 = 4000   # deg K
P0 = 1e9    # start of pressure interval, Pa
P1 = 100e9  # End of pressure interval, Pa
dT = 50      # Temperature step
dP = 1e9    # Pressure step
T = arange(T0, T1, dT)
P = arange(P0, P1, dP)
nT = len(T)-1
nP = len(P) - 1
m = 0.0403044 # molar mass of MgO, kg/mol

# R = 8.314 # gas constant, j/mol
# Pr = 1e5 # reference pressure, 1 bar
# Tr = 298.15 # reference plot_surfaceace temperature, K
# #
# # empirical parameters
# Hr      = -6.015e5    # J empirical constant
# Vr      = 1.12228e-5  # J/Pa
# phi     = 3.01795e10 # Pa
# c1      = 1.96612
# c2      = 4.12756
# c3      = 0.53690
# deltaH1 = 2966.88 # J
# deltaH2 = 56212.69 # J
# deltaH3 = 27787.19 # J
# deltaV1 = deltaV2 = 3.52971e-8 # J/Pa
# deltaV3 = 1.9849568e-6 # J/Pa


def V(press, temp):
    return (Gibbs_m(press + dP, temp) - Gibbs_m(press, temp)) / dP

def rho(press, temp):
    return m / V(press, temp)

def alpha(press, temp):
    return (-1/rho(press, temp)) * (rho(press, temp + dT) - rho(press, temp)) / dT

def beta(press, temp):
    return 1/rho(press,temp) * (rho(press + dP, temp) - rho(press, temp)) / dP



# visualise
fig, axes = plt.subplots(3,3, figsize=(20,12))
volume = density = expansivity = compressibility = zeros((nT, nP))
for i in range(nT):
    for j in range(nP):
        volume[i,j] = V(press=P[j], temp=T[i])
        density[i,j] = rho(press=P[j], temp=T[i])
        expansivity[i,j] = alpha(press=P[j], temp=T[i])
        compressibility[i,j] = beta(press=P[j], temp=T[i])
    print(T[i]) if i % 10 == 0 else None

PP, TT = meshgrid(P[0:-1]/1e9, T[0:-1])

# pcolormeshes
for i in range(3):
    if i == 0:
        pc = axes[i][0].pcolormesh(TT, PP, density)
        fig.colorbar(pc, ax=axes[i][0])
        cf = axes[i][1].contourf(TT, PP, density)
        ax = fig.add_subplot(3,3,i*3+3, projection='3d')
        s = ax.plot_surface(TT, PP, density)
    if i == 1:
        pc = axes[i][0].pcolormesh(TT, PP, expansivity)
        fig.colorbar(pc, ax=axes[i][0])
        cf = axes[i][1].contourf(TT, PP, expansivity)
        ax = fig.add_subplot(3,3,i*3+3, projection='3d')

        s = ax.plot_surface(TT, PP, expansivity)
    if i == 2:
        pc = axes[i][0].pcolormesh(TT, PP, compressibility)
        fig.colorbar(pc, ax=axes[i][0])
        cf = axes[i][1].contourf(TT, PP, compressibility)
        ax = fig.add_subplot(3,3,i*3+3, projection='3d')

        s = ax.plot_surface(TT, PP, compressibility)
        axes[i][0].set_xlabel("Temperature [K]")
        axes[i][0].set_ylabel("Pressure [GPa]")
    plt.pause(0.001)


plt.show()
