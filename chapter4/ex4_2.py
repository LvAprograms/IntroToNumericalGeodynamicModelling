from matplotlib import pyplot as plt
from math import sqrt
from numpy import zeros, array, meshgrid, linspace, sin, cos, pi
# from chapter1.ex1 import make_grid, calculate_properties

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
    eii = sqrt(0.5 * (exx ** 2 + eyy ** 2) + exy ** 2)
    return exx, eyy, exy, eii

    
W = 1E6    # width of box in m
H = 1.5E7  # height of box in m
nx = 31    # nodes in x direction
ny = 31    # nodes in y direction
x,y = make_grid(nx, ny, W, H)
xx, yy = meshgrid(x, y)

vx0 = 1e-9 * W/2/H # scaling factor for horizonta velocity, m/s
vy0 = 1e-9         # scaling factor for vertical velocity, m/s

# calculate vx, vy, dvx/dx, dvy/dy, div(v) on the grid
vx, vy, dvxdx, dvydy, div_v, dvxdy, dvydx = calculate_properties(vx0, vy0, xx, yy, W, H)
exx, eyy, exy, eii = calc_dev_strainrates(vx, vy, x, y, dvxdx, dvydy, dvxdy, dvydx)




# plot properties
fig, axes = plt.subplots(2, 3, figsize=(20,10))
plt.set_cmap('hot')
properties = [vx, vy, dvxdx, dvydy, div_v, eii]
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