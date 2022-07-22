from numpy import zeros, arange, meshgrid, mean, reshape, array
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from gibbsfree import Gibbs_m


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


def V(press, temp):
    a = Gibbs_m(press + dP, temp)
    b = Gibbs_m(press, temp)
    return (a - b) / dP


def rho(press, temp):
    return m / V(press, temp)


def alpha(press, temp):
    return (-1/rho(press,temp)) * ((rho(press, temp + dT) - rho(press, temp)) / dT)


def beta(press, temp):
    return 1/rho(press,temp) * ((rho(press + dP, temp) - rho(press, temp)) / dP)


# visualise
fig = plt.figure()
volume = zeros((nT, nP))
density = zeros((nT, nP))
expansivity = zeros((nT, nP))
compressibility = zeros((nT, nP))
for i in range(nT):
    for j in range(nP):
        volume[i,j] = V(press=P[j], temp=T[i])
        density[i,j] = m / volume[i,j]
        expansivity[i,j] = alpha(press=P[j], temp=T[i])
        compressibility[i,j] = beta(press=P[j], temp=T[i])
    print([mean(volume), mean(density), mean(expansivity), mean(compressibility)]) if i % 10 == 0 else None

PP, TT = meshgrid(P[0:-1]/1e9, T[0:-1])
axlist = []
# pcolormeshes
for i in range(3):
    if i == 0:
        ax = fig.add_subplot(3,3,i*3+1)
        axlist.append(ax)
        pc = ax.pcolormesh(TT, PP, density, cmap='hot')
        fig.colorbar(pc, ax=ax)
        ax = fig.add_subplot(3,3,i*3+2)
        axlist.append(ax)
        cf = ax.contourf(TT, PP, density, cmap='hot')
        ax = fig.add_subplot(3,3,i*3+3, projection='3d')
        axlist.append(ax)
        s = ax.plot_surface(TT, PP, density, cmap='hot')
    if i == 1:
        ax = fig.add_subplot(3, 3, i * 3 + 1)
        pc = ax.pcolormesh(TT, PP, expansivity, cmap='hot')
        fig.colorbar(pc, ax=ax)
        axlist.append(ax)
        ax = fig.add_subplot(3, 3, i * 3 + 2)
        axlist.append(ax)
        cf = ax.contourf(TT, PP, expansivity, cmap='hot')
        ax = fig.add_subplot(3, 3, i * 3 + 3, projection='3d')
        axlist.append(ax)
        s = ax.plot_surface(TT, PP, expansivity, cmap='hot')
    if i == 2:
        ax = fig.add_subplot(3, 3, i * 3 + 1)
        axlist.append(ax)
        pc = ax.pcolormesh(TT, PP, compressibility, cmap='hot')
        fig.colorbar(pc, ax=ax)
        ax = fig.add_subplot(3, 3, i * 3 + 2)
        axlist.append(ax)
        cf = ax.contourf(TT, PP, compressibility, cmap='hot')
        ax = fig.add_subplot(3, 3, i * 3 + 3, projection='3d')
        axlist.append(ax)
        s = ax.plot_surface(TT, PP, compressibility, cmap='hot')

    plt.pause(0.001)
titles = [r'$\rho$', r'$\alpha$', r'$\beta$']
plottype = [' pcolor', ' contourf', ' surf']
for i in range(3):
    for j in range(3):
        axlist[i * 3 + j].set_xlabel("Temperature [K]")
        axlist[i * 3 + j].set_ylabel("Pressure [GPa]")
        axlist[i * 3 + j].set_title(titles[i] + plottype[j])
plt.show()

## Exercise 2.3
# load from files
nx = 0
ny = 0
T0 = 0
P0 = 0
T_step = 0
P_step = 0
rho_pyrolite = []
with open('m895_ro', 'r') as f:
    for i, line in enumerate(f.readlines()):
        if i > 0:
            l = line.strip('\n').split()
            if i == 1:
                nx = int(l[0])
                ny = int(l[1])
                T0 = float(l[2])
                P0 = float(l[3])
                T_step = float(l[4])
                P_step = float(l[5])
            elif i > 2:
                rho_pyrolite.append(float(l[0]))

rho_pyrolite = reshape(array(rho_pyrolite), (nx, ny))
T = [T0 + (T_step * i) for i in range(nx)]
P = [P0 + (P_step * i) for i in range(ny)]
TT, PP = meshgrid(T, P)

rho_morb = []
with open('morn_ro', 'r') as f:
    for i, line in enumerate(f.readlines()):
        if i > 0:
            l = line.strip('\n').split()
            if i == 1:
                nx = int(l[0])
                ny = int(l[1])
                T0 = float(l[2])
                P0 = float(l[3])
                T_step = float(l[4])
                P_step = float(l[5])
            elif i > 2:
                rho_morb.append(float(l[0]))

rho_morb = reshape(array(rho_morb), (nx, ny))
## visualizing
fig1, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(24, 8))
pc1 = ax1.pcolormesh(TT, PP, rho_pyrolite, cmap='hot')
fig1.colorbar(pc1, ax=ax1)
pc2 = ax2.pcolormesh(TT, PP, rho_morb, cmap='hot')
fig1.colorbar(pc2, ax=ax2)
pc3 = ax3.pcolormesh(TT, PP, rho_pyrolite - rho_morb, cmap='hot')
fig1.colorbar(pc3, ax=ax3)


plt.show()



