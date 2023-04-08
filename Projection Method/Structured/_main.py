import os
import sys
import math as m
import numpy as np
import matplotlib.pyplot as plt

# domain parameters
L = 3.0  # channel length [m]
H = 1.0  # channel height [m]
nx = 150  # cells in x-direction (even/odd)
ny = 100 + 1  # cells in y-direction (odd)
dx = L/nx  # cell dimension in x-direction [m]
dy = H/ny  # cell dimension in y-direction [m]

a = 25 + 1  # cell position in x-direction of square cylinder centre
b = int(ny/2) + 1  # cell position in y-direction of square cylinder centre
cx = 5  # cells in x-direction from square cylinder centre
cy = 10  # cells in y-direction from square cylinder centre

# flow parameters
U = 5.0  # inlet velocity [m/s]
nu = 0.01  # kinematic viscosity [m2/s]
Re = U*(cy+cy+1)/ny*H/nu  # Reynolds number

t = 10.0  # total time [s]
dl = min(dx, dy)  # minimum cell dimension [m]
dt = min(0.1*dl/U, 0.1*dl*dl/nu)  # time step [s]
nt = int(t/dt)  # time iterations

print('Re =', Re)
print('dt =', dt)

# output parameters
data = 'data'  # data folder name
interval = 100  # save interval
os.mkdir(data)  # data directory

figure = 'figure'  # figure folder name
limit = 15.0  # vorticity limit [/s]
os.mkdir(figure)  # figure directory

###################
## program start ##
###################

# flow variables
u = np.zeros((nx + 2, ny + 2))  # cell u-velocity
v = np.zeros((nx + 2, ny + 2))  # cell v-velocity
p = np.zeros((nx + 2, ny + 2))  # cell pressure

vN = np.zeros((nx + 2, ny + 2))  # north face mass flow rate
vS = np.zeros((nx + 2, ny + 2))  # south face mass flow rate
uE = np.zeros((nx + 2, ny + 2))  # east face mass flow rate
uW = np.zeros((nx + 2, ny + 2))  # west face mass flow rate

# u-velocity control volume fluxes
FudN = np.zeros((nx + 2, ny + 2))  # north face diffusion flux
FudS = np.zeros((nx + 2, ny + 2))  # south face diffusion flux
FudE = np.zeros((nx + 2, ny + 2))  # east face diffusion flux
FudW = np.zeros((nx + 2, ny + 2))  # west face diffusion flux

FucN = np.zeros((nx + 2, ny + 2))  # north face convection flux
FucS = np.zeros((nx + 2, ny + 2))  # south face convection flux
FucE = np.zeros((nx + 2, ny + 2))  # east face convection flux
FucW = np.zeros((nx + 2, ny + 2))  # west face convection flux

# v-velocity control volume fluxes
FvdN = np.zeros((nx + 2, ny + 2))  # north face diffusion flux
FvdS = np.zeros((nx + 2, ny + 2))  # south face diffusion flux
FvdE = np.zeros((nx + 2, ny + 2))  # east face diffusion flux
FvdW = np.zeros((nx + 2, ny + 2))  # west face diffusion flux

FvcN = np.zeros((nx + 2, ny + 2))  # north face convection flux
FvcS = np.zeros((nx + 2, ny + 2))  # south face convection flux
FvcE = np.zeros((nx + 2, ny + 2))  # east face convection flux
FvcW = np.zeros((nx + 2, ny + 2))  # west face convection flux

# pressure poisson equation coefficients
aN = (1/dy/dy)*np.ones((nx + 2, ny + 2))  # north cell coefficient
aS = (1/dy/dy)*np.ones((nx + 2, ny + 2))  # south cell coefficient
aE = (1/dx/dx)*np.ones((nx + 2, ny + 2))  # east cell coefficient
aW = (1/dx/dx)*np.ones((nx + 2, ny + 2))  # west cell coefficient
aP = np.zeros((nx + 2, ny + 2))  # current cell coefficient
bP = np.zeros((nx + 2, ny + 2))  # current cell source term

aN[1:-1, -2] = 0.0  # top wall
aS[1:-1, 1] = 0.0  # bottom wall
aE[-2, 1:-1] = 0.0  # outlet
aW[1, 1:-1] = 0.0  # inlet

aN[(a - cx):(a + cx + 1), b - cy - 1] = 0.0  # cylinder south
aS[(a - cx):(a + cx + 1), b + cy + 1] = 0.0  # cylinder north
aE[a - cx - 1, (b - cy):(b + cy + 1)] = 0.0  # cylinder west
aW[a + cx + 1, (b - cy):(b + cy + 1)] = 0.0  # cylinder east

aP[1:-1, 1:-1] = -(aN[1:-1, 1:-1] + aS[1:-1, 1:-1] + aE[1:-1, 1:-1] + aW[1:-1, 1:-1])

# boundary conditions
u[1, 1:-1] = U  # inlet

err = np.zeros(nt)
for n in range(nt):
    
    ###############################
    ## u-velocity control volume ##
    ###############################

    # mass flow rate
    vN[2:-1, 1:-1] = 0.5*(v[1:-2, 2:] + v[2:-1, 2:])
    vS[2:-1, 1:-1] = 0.5*(v[1:-2, 1:-1] + v[2:-1, 1:-1])
    uE[2:-1, 1:-1] = 0.5*(u[2:-1, 1:-1] + u[3:, 1:-1])
    uW[2:-1, 1:-1] = 0.5*(u[2:-1, 1:-1] + u[1:-2, 1:-1])

    # north flux
    # top wall
    FudN[2:-1, -2] = -nu*(9*u[2:-1, -2] - u[2:-1, -3])/(3*dy)*dx
    FucN[2:-1, -2] = 0.0
    # inner cells
    FudN[2:-1, 1:-2] = nu*(u[2:-1, 2:-1] - u[2:-1, 1:-2])/dy*dx
    FucN[2:-1, 1:-2] = (np.maximum(vN[2:-1, 1:-2], 0)*u[2:-1, 1:-2] - np.maximum(-vN[2:-1, 1:-2], 0)*u[2:-1, 2:-1])*dx
    # square cylinder south
    FudN[(a - cx):(a + cx + 2), b - cy - 1] = -nu*(9*u[(a - cx):(a + cx + 2), b - cy - 1] - u[(a - cx):(a + cx + 2), b - cy - 2])/(3*dy)*dx
    FucN[(a - cx):(a + cx + 2), b - cy - 1] = 0.0

    # south flux
    # bottom wall
    FudS[2:-1, 1] = nu*(9*u[2:-1, 1] - u[2:-1, 2])/(3*dy)*dx
    FucS[2:-1, 1] = 0.0
    # inner cells
    FudS[2:-1, 2:-1] = nu*(u[2:-1, 2:-1] - u[2:-1, 1:-2])/dy*dx
    FucS[2:-1, 2:-1] = (np.maximum(vS[2:-1, 2:-1], 0)*u[2:-1, 1:-2] - np.maximum(-vS[2:-1, 2:-1], 0)*u[2:-1, 2:-1])*dx
    # square cylinder north
    FudS[(a - cx):(a + cx + 2), b + cy + 1] = nu*(9*u[(a - cx):(a + cx + 2), b + cy + 1] - u[(a - cx):(a + cx + 2), b + cy + 2])/(3*dy)*dx
    FucS[(a - cx):(a + cx + 2), b + cy + 1] = 0.0
    
    # east flux
    # outlet
    FudE[-2, 1:-1] = 0.0
    FucE[-2, 1:-1] = (u[-1, 1:-1]*u[-1, 1:-1])*dy
    # inner cells
    FudE[2:-2, 1:-1] = nu*(u[3:-1, 1:-1] - u[2:-2, 1:-1])/dx*dy
    FucE[2:-2, 1:-1] = (np.maximum(uE[2:-2, 1:-1], 0)*u[2:-2, 1:-1] - np.maximum(-uE[2:-2, 1:-1], 0)*u[3:-1, 1:-1])*dy

    # west flux
    # inlet
    FudW[2, 1:-1] = nu*(u[2, 1:-1] - U)/dx*dy
    FucW[2, 1:-1] = U*U*dy
    # inner cells
    FudW[3:-1, 1:-1] = nu*(u[3:-1, 1:-1] - u[2:-2, 1:-1])/dx*dy
    FucW[3:-1, 1:-1] = (np.maximum(uW[3:-1, 1:-1], 0)*u[2:-2, 1:-1] - np.maximum(-uW[3:-1, 1:-1], 0)*u[3:-1, 1:-1])*dy
    
    ###############################
    ## v-velocity control volume ##
    ###############################

    # mass flow rate
    vN[1:-1, 2:-1] = 0.5*(v[1:-1, 2:-1] + v[1:-1, 3:])
    vS[1:-1, 2:-1] = 0.5*(v[1:-1, 2:-1] + v[1:-1, 1:-2])
    uE[1:-1, 2:-1] = 0.5*(u[2:, 1:-2] + u[2:, 2:-1])
    uW[1:-1, 2:-1] = 0.5*(u[1:-1, 2:-1] + u[1:-1, 1:-2])

    # north flux
    # top wall
    FvdN[1:-1, -2] = -nu*v[1:-1, -2]/dy*dx
    FvcN[1:-1, -2] = (np.maximum(vN[1:-1, -2], 0)*v[1:-1, -2] - np.maximum(-vN[1:-1, -2], 0)*v[1:-1, -1])*dx
    # inner cells
    FvdN[1:-1, 2:-2] = nu*(v[1:-1, 3:-1] - v[1:-1, 2:-2])/dy*dx
    FvcN[1:-1, 2:-2] = (np.maximum(vN[1:-1, 2:-2], 0)*v[1:-1, 2:-2] - np.maximum(-vN[1:-1, 2:-2], 0)*v[1:-1, 3:-1])*dx

    # south flux
    # bottom wall
    FvdS[1:-1, 2] = nu*v[1:-1, 2]/dy*dx
    FvcS[1:-1, 2] = (np.maximum(vS[1:-1, 2], 0)*v[1:-1, 1] - np.maximum(-vS[1:-1, 2], 0)*v[1:-1, 2])*dx
    # inner cells
    FvdS[1:-1, 3:-1] = nu*(v[1:-1, 3:-1] - v[1:-1, 2:-2])/dy*dx
    FvcS[1:-1, 3:-1] = (np.maximum(vS[1:-1, 3:-1], 0)*v[1:-1, 2:-2] - np.maximum(-vS[1:-1, 3:-1], 0)*v[1:-1, 3:-1])*dx

    # east flux
    # outlet
    FvdE[-2, 2:-1] = 0.0
    FvcE[-2, 2:-1] = (np.maximum(uE[-2, 2:-1], 0)*v[-2, 2:-1] - np.maximum(-uE[-2, 2:-1], 0)*v[-1, 2:-1])*dy
    # inner cells
    FvdE[1:-2, 2:-1] = nu*(v[2:-1, 2:-1] - v[1:-2, 2:-1])/dx*dy
    FvcE[1:-2, 2:-1] = (np.maximum(uE[1:-2, 2:-1], 0)*v[1:-2, 2:-1] - np.maximum(-uE[1:-2, 2:-1], 0)*v[2:-1, 2:-1])*dy
    # square cylinder west
    FvdE[a - cx - 1, (b - cy):(b + cy + 2)] = -nu*(9*v[a - cx - 1, (b - cy):(b + cy + 2)] - v[a - cx - 2, (b - cy):(b + cy + 2)])/(3*dx)*dy
    FvcE[a - cx - 1, (b - cy):(b + cy + 2)] = 0.0

    # west flux
    # inlet
    FvdW[1, 2:-1] = nu*v[1, 2:-1]/dx*dy
    FvcW[1, 2:-1] = 0.0
    # inner cells
    FvdW[2:-1, 2:-1] = nu*(v[2:-1, 2:-1] - v[1:-2, 2:-1])/dx*dy
    FvcW[2:-1, 2:-1] = (np.maximum(uW[2:-1, 2:-1], 0)*v[1:-2, 2:-1] - np.maximum(-uW[2:-1, 2:-1], 0)*v[2:-1, 2:-1])*dy
    # square cylinder east
    FvdW[a + cx + 1, (b - cy):(b + cy + 2)] = nu*(9*v[a + cx + 1, (b - cy):(b + cy + 2)] - v[a + cx + 2, (b - cy):(b + cy + 2)])/(3*dx)*dy
    FvcW[a + cx + 1, (b - cy):(b + cy + 2)] = 0.0

    #########################
    ## velocity projection ##
    #########################

    # intermediate velocity
    Fuc = FucN - FucS + FucE - FucW
    Fud = FudN - FudS + FudE - FudW
    Fvc = FvcN - FvcS + FvcE - FvcW
    Fvd = FvdN - FvdS + FvdE - FvdW
    u = u + dt/(dx*dy)*(-Fuc + Fud)
    v = v + dt/(dx*dy)*(-Fvc + Fvd)

    # boundary conditions
    u[-1, 1:-1] = u[-2, 1:-1]  # outlet
    v[-1, 2:-1] = v[-2, 2:-1]  # outlet
    u[(a - cx):(a + cx + 2), (b - cy):(b + cy + 1)] = 0.0  # square cylinder
    v[(a - cx):(a + cx + 1), (b - cy):(b + cy + 2)] = 0.0  # square cylinder

    #########################
    ## pressure correction ##
    #########################

    # pressure poisson
    bP[1:-1, 1:-1] = (1/dt)*((u[2:, 1:-1] - u[1:-1, 1:-1])/dx + (v[1:-1, 2:] - v[1:-1, 1:-1])/dy)
    for k in range(10):  # Gauss-Seidel
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                p[i, j] = (1/aP[i, j])*(-aN[i, j]*p[i, j + 1] - aS[i, j]*p[i, j - 1] - aE[i, j]*p[i + 1, j] - aW[i, j]*p[i - 1, j] + bP[i, j])
    
    # pressure correction
    u[2:-1, 1:-1] = u[2:-1, 1:-1] - (dt/dx)*(p[2:-1, 1:-1] - p[1:-2, 1:-1])
    v[1:-1, 2:-1] = v[1:-1, 2:-1] - (dt/dy)*(p[1:-1, 2:-1] - p[1:-1, 1:-2])

    # boundary conditions
    # u[-1, 1:-1] = u[-2, 1:-1]  # outlet
    u[-1, 1:-1] = u[-2, 1:-1] - (dx/dy)*(v[-2, 2:] - v[-2, 1:-1])  # outlet
    v[-1, 2:-1] = v[-2, 2:-1]  # outlet
    u[(a - cx):(a + cx + 2), (b - cy):(b + cy + 1)] = 0.0  # square cylinder
    v[(a - cx):(a + cx + 1), (b - cy):(b + cy + 2)] = 0.0  # square cylinder

    ######################
    ## check continuity ##
    ######################

    err[n] = (sum(u[-1, 1:-1]) - sum(u[1, 1:-1]))/sum(u[1, 1:-1])*100
    print('Iteration = %d of %d (%.5f%%), Global Mass Imbalance [%%] = %.5f' % (n, nt - 1, n/(nt - 1)*100, err[n]))

    ###############
    ## save data ##
    ###############

    if n % interval == 0:
        np.save(os.path.join(data, 'u{:d}'.format(n)), u)
        np.save(os.path.join(data, 'v{:d}'.format(n)), v)
        np.save(os.path.join(data, 'p{:d}'.format(n)), p)

###############
## plot data ##
###############

# plot global mass imbalance
plt.plot(np.linspace(0, t, nt), err)
plt.title('Global Mass Imbalance', fontweight='bold')
plt.xlabel('Time [s]')
plt.ylabel('Imbalance [%]')
plt.grid('on')
plt.show()

# plot vorticity
plt.figure(figsize=(15, 5))
x = np.linspace(0, L, nx)
y = np.linspace(0, H, ny)
xx, yy = np.meshgrid(x, y)

for n in range(0, nt, interval):
    u = np.load(os.path.join(sys.path[0], data, 'u{:d}.npy'.format(n)))
    v = np.load(os.path.join(sys.path[0], data, 'v{:d}.npy'.format(n)))
    dvdx = (v[2:, 2:] + v[2:, 1:-1] - v[:-2, 2:] - v[:-2, 1:-1])/(4*dx)
    dudy = (u[1:-1, 2:] + u[2:, 2:] - u[1:-1, 1:-1] - u[2:, 1:-1])/(4*dy)
    w = dvdx - dudy

    cmap = plt.get_cmap('viridis')
    levels = np.linspace(-limit, limit, int(4*limit) + 1)
    plt.contourf(xx, yy, np.transpose(w), cmap=cmap, levels=levels, extend='both')
    plt.colorbar(label='Vorticity [/s]')

    sx = np.array([a - cx - 1, a - cx - 1, a + cx - 1, a + cx - 1, a - cx - 1])/nx*L
    sy = np.array([b - cy - 1, b + cy - 1, b + cy - 1, b - cy - 1, b - cy - 1])/ny*H
    plt.fill(sx, sy, zorder=10)

    plt.title('Re = %d, t = %.3f s' % (Re, n*dt), fontweight='bold')
    plt.xlabel('Channel Length [m]')
    plt.ylabel('Channel Height [m]')
    plt.axis('equal')
    plt.savefig(os.path.join(sys.path[0], figure, 'figure{:d}.png'.format(n)))
    plt.clf()
    print(n)
