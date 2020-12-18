from matplotlib import pyplot as plt
from numpy import arange

def one_D_steadyStateGeotherm(k=2, Hr=1E-6, zbase=35e3, T0=300, T1=700, figure=None):
    """
    Plots a 1D steady state continental geotherm. the `figure` param can be specified if a parameter study is the goal. Otherwise, only one plot is created.
    """
    z = arange(0, zbase+1e3, 1e3)
    c1 = T0
    c0 = (T1- T0) / zbase + zbase * (Hr/(2*k))

    T = -Hr / (2*k) * z**2 + c0 * z + c1
    if figure is None:
        fig, ax = plt.subplots()
        plt.gca().invert_yaxis()

    else:
        fig = figure[0]
        ax = figure[1]

    ax.plot(T, z, label="k={},Hr={:.2f},T0={},T1={}".format(k,Hr*1e6,T0, T1))
    fig.legend()
    plt.grid(b=True)
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("Depth [m]")
    plt.pause(0.001)

if __name__ == "__main__":
    fig, ax = plt.subplots()
    Hr=[0, 1, 2, 5]
    for v in Hr:
        one_D_steadyStateGeotherm(Hr=v*1e-6, figure=(fig, ax))        
    plt.gca().invert_yaxis()

    plt.show()