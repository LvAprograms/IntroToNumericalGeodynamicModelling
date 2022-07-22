from matplotlib import pyplot as plt
from math import sqrt, log10
from numpy import zeros, array, meshgrid, linspace, sin, cos, pi, reshape
# from chapter1.ex1 import make_grid, calculate_properties

# FUNCTION DEFINITIONS -----------------------------------------------------------------------
# redefined from ex 1.1 because of VScode :( 
def make_grid(nx=31, ny=31, width=1e6, height=15e6):
    x = linspace(0, width, num=nx)
    y = linspace(0, height, num=ny)
    return x, y

def calculate_properties(vx0, vy0, x, y, width, height):
    vx = -vx0 * sin(2 * pi * x / width) * cos(pi * y / height)
    vy = vy0 * cos(2 * pi * x / width) * sin(pi * y / height)
    dvxdx = -vx0 * (2 * pi / width) * cos(2 * pi * x / width) * cos(pi * y / height)
    dvydy = vy0 * (pi / height) * cos(2 * pi * x / width) * cos(pi * y / height)
    dvxdy = pi * vx0 / height * sin(2 * pi * x / width) * sin(pi * y / height)
    dvydx = -2 * pi / width * vy0 * sin(2 * pi * x / width)  * sin(pi * y / height)
    div_v = dvxdx + dvydy
    return vx, vy, dvxdx, dvydy, div_v, dvxdy, dvydx

def calc_dev_strainrates(vx, vy, x, y, dvxdx, dvydy, dvxdy, dvydx):
    trace=dvxdx + dvydy
    exy = 0.5 * (dvxdy + dvydx) - 0.5 * trace
    exx = dvxdx - 0.5 * trace
    eyy = dvydy - 0.5 * trace
    eii = []
    for i in range(nx):
        for j in range(ny):
            eii.append(sqrt(0.5 * (exx[i][j] ** 2 + eyy[i][j] ** 2) + exy[i][j] ** 2))
    return exx, eyy, exy, reshape(array(eii), (nx, ny))

if __name__ == "__main__":
    # Model variable definitions -------------------------------------------------------------
    W = 1E6    # width of box in m
    H = 1.5E7  # height of box in m
    nx = 31    # nodes in x direction
    ny = 31    # nodes in y direction
    x,y = make_grid(nx, ny, W, H)
    xx, yy = meshgrid(x, y)

    vx0 = 1e-9 * W/2/H # scaling factor for horizonta velocity, m/s
    vy0 = 1e-9         # scaling factor for vertical velocity, m/s

    # Calculating the needed variables -------------------------------------------------
    vx, vy, dvxdx, dvydy, div_v, dvxdy, dvydx = calculate_properties(vx0, vy0, xx, yy, W, H)
    exx, eyy, exy, eii = calc_dev_strainrates(vx, vy, x, y, dvxdx, dvydy, dvxdy, dvydx)


    # Plot the variables with LaTeX formattied labels--------------------------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(20,10))
    plt.set_cmap('hot')
    properties = [vx, vy, exx, eyy, exy, eii]
    labels = ["Vx", "Vy", r'$\dot{\epsilon^{\prime}}_{xx}$',  r'$\dot{\epsilon^{\prime}}_{yy}$', r'$\dot{\epsilon^{\prime}}_{xy}$', r'$\dot{\epsilon\prime}_{II}$']
    for i, p in enumerate(properties):
        ax = axes.flatten()[i]
        pc = ax.pcolormesh(xx, yy, p, shading='auto' )
        ax.quiver(xx,yy, -vx, -vy) # overlay everything with velocity arrows
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m'")
        ax.set_title(labels[i])
        fig.colorbar(pc, ax=ax)

    plt.savefig("chapter4/ex4_2.png")