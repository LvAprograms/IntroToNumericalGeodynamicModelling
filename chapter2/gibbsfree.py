from numpy import log, exp

def Gibbs_m(P, T):
    """
    Returns gibbs free energy as a function of temperature and pressure
    """
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
    deltaH2 = 5621.69 # J
    deltaH3 = 27787.19 # J
    deltaV1 = deltaV2 = 3.52971e-8 # J/Pa
    deltaV3 = 1.9849568e-6 # J/Pa
    c = [c1, c2, c3]
    dH = [deltaH1, deltaH2, deltaH3]
    dV = [deltaV1, deltaV2, deltaV3]

    Psi = 0
    sumpart = 0
    for i in range(3):
        Psi = (5/4) * (Pr + phi) **(1/5) * ((P + phi)**(4/5) - (Pr + phi)**(4/5))
        e = exp(-((dH[i] + dV[i] * Psi))/(R * T)) 
        eo = exp(-dH[i] / (R * Tr))
        sumpart += c[i] * (R * T * log(1 - e) - (dH[i] * eo / (1 - eo)))
    return Hr + Vr * Psi + sumpart
