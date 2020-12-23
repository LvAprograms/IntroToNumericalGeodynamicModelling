from numpy import array, arange
from matplotlib import pyplot as plt

## Exercise 2.2 (2.1 is analytical)
T0 = 100    # deg K
T1 = 4000   # deg K
P0 = 1e9    # start of pressure interval, Pa
P1 = 100e9  # End of pressure interval, Pa
T = arange(T0, T1, 1)
P = arange(P0, P1, 1E9)

R = 8.314 # gas constant, j/mol
Pr = 1e5 # reference pressure, 1 bar
Tr = 298.15 # reference surface temperature, K

# empirical parameters
Hr      = -6.015e5    # J empirical constant
Vr      = 1.12228e-5  # J/Pa
phi     = 3.01795e10 # Pa
c1      = 1.96612 
c2      = 4.12756
c3      = 0.53690
deltaH1 = 2966.88 # J
deltaH2 = 56212.69 # J
deltaH3 = 27787.19 # J
deltaV1 = deltaV2 = 3.52971e-8 # J/Pa
deltaV3 = 1.9849568e-6 # J/Pa


