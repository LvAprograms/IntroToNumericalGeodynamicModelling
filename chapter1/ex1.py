from numpy import linspace, sin, cos, pi, array, meshgrid
from matplotlib import pyplot as plt
# exercise 1.1

def v_mantle(x, y, z):
    v = []
    vx = 1E-10 + 1e-13 * (x+y+z)
    vy = 1e-10 - 1e-13 * (x - 2*y - 3*z)
    vz = 1e-10 - 1e-13 * (x + y + 2*z)
    v = [vx, vy, vz]
    return v

def rho_mantle(x, y, z):
    return 3300 + 0.001*(x - 2*y + z)

def ex1_1():
    x, y, z = 1000, 1000, 1000
    v = v_mantle(x,y,z)
    rho = rho_mantle(x, y, z)
    # div_v = dvx/dx + dvy/dz + dvz/dz = 1e-13 + 2e-13 -2e-13 = 1e-13
    # drho/dt = 0
    # Drho/Dt = drho/dt + v * grad(rho)
    DrhoDt = v[0] * 0.001 + v[1] * -0.002 + v[2] * 0.001
    print("rho = {}, div(v) = 1e-13, drho/dt = 0 and Drho/Dt = {}".format(rho, DrhoDt))

# exercise 1.2
# write a matlab ----> python :P code for computing and visualising a 2D velocity field and its divergence. Model design:
# an area of mantle (1000x1500km2) is convecting with one central upwelling in the middle of the model box and two downwellings at
# the sides. 
# vx = -vx0 * sin(2 * pi * x/W) * cos(pi * y/H)
# vy = vy0 * cos(2* pi *x/W) * sin(pi * y/H)

def make_grid(nx=31, ny=31, width=1e6, height=15e6):
    x = linspace(0, width, num=nx)
    y = linspace(0, height, num=ny)
    return x, y

def calculate_properties(vx0, vy0, x, y, width, height):
    vx = -vx0 * sin(2 * pi * x / width) * cos(pi * y / height)
    vy = vy0 * cos(2 * pi * x / width) * sin(pi * y / height)
    dvxdx = -vx0 * (2 * pi / width) * cos(2 * pi * x / width) * cos(pi * y / height)
    dvydy = vy0 * (pi / height) * cos(2 * pi * x / width) * cos(pi * y / height)
    div_v = dvxdx + dvydy
    return vx, vy, dvxdx, dvydy, div_v

def ex1_2():
    W = 1E6    # width of box in m
    H = 1.5E7  # height of box in m
    nx = 31    # nodes in x direction
    ny = 31    # nodes in y direction
    x,y = make_grid(nx, ny, W, H)
    xx, yy = meshgrid(x, y)

    vx0 = 1e-9 * W/2/H # scaling factor for horizonta velocity, m/s
    vy0 = 1e-9         # scaling factor for vertical velocity, m/s
    
    # calculate vx, vy, dvx/dx, dvy/dy, div(v) on the grid
    vx, vy, dvxdx, dvydy, div_v = calculate_properties(vx0, vy0, xx, yy, W, H)

    # plot properties
    fig, axes = plt.subplots(2, 3, figsize=(20,10))
    plt.set_cmap('hot')
    properties = [vx, vy, dvxdx, dvydy, div_v]
    labels = ["Vx", "Vy", "dVx/dx", "dVy/dy", "div(v)", "arrows?"]
    for i, p in enumerate(properties):
        ax = axes.flatten()[i]
        pc = ax.pcolormesh(xx, yy, p, shading='interp' )
        ax.quiver(xx,yy, -vx, -vy) # overlay everything with velocity arrows
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m'")
        ax.set_title(labels[i])
        fig.colorbar(pc, ax=ax)
    plt.show()

if __name__ == "__main__":
    ex1_1()
    ex1_2()