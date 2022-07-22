from matplotlib import pyplot as plt
from numpy import zeros, meshgrid



# We have to read two files with the same format, so I define a function that reads and processes the data. 
def read_densities(file, TP=False):
    # initialise variables
    dens = zeros((1,1)) 
    nT, nP = 0, 0
    Tstart, Pstart, Tstep, Pstep = 0.0, 0.0, 0.0, 0.0
    T, P = zeros((nP, nT)), zeros((nP, nT))
    # read the file and fill density matrix
    with open(file) as f:
        T, P = [], []
        for i, line in enumerate(f.readlines()):
            l = line.split()
            if i == 1:
                nT = int(l[0])
                nP = int(l[1])
                dens = zeros((nT, nP))
                Tstart = float(l[2])
                Pstart = float(l[3])
                Tstep = float(l[4])
            elif i == 2: 
                Pstep = float(l[0])
            elif i > 3:
                dens[int((i-4)/nP), (i-4)%nT] = float(l[0]) # data starts at i = 4, each row in the matrix is nT long, so every nTth entry the row will increase. 
    print("Tstart: {} K, Pstart: {} bar, Tstep: {} K, Pstep: {} bar".format(Tstart, Pstart, Tstep, Pstep))

    if TP: # only calculate Temperature and Density vectors if not already done before
        T = [Tstart + _ * Tstep for _ in range(nT)]
        P = [Pstart + _ * Pstep for _ in range(nP)] 
        return meshgrid(T, P), dens
    else:
        return dens 
                
            
[T, P], pyrolite = read_densities("chapter2/m895_ro", TP=True)
morb = read_densities("chapter2/morn_ro")
rhodiff = morb - pyrolite
plt.figure()
plt.contourf(T, P, rhodiff, cmap='hot')
plt.xlabel("Temperature [K]")
plt.ylabel("Pressure [bar]")
c = plt.colorbar()
c.set_label(r'$\Delta\rho$ [kg/m3]')
plt.show()