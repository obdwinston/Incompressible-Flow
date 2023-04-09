import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import moviepy.editor as mpy

folder = 'data/'
variable = 'speed'
interval = 50

U = 5.  # inlet x-velocity [m/s]
V = 0.  # inlet y-velocity [m/s]
nu = .01  # kinematic viscosity [m2/s]
t = 5.  # time [s]
dt = 1e-4  # time step [s]

L = 3.  # domain length [m]
H = 1.  # domain height [m]
D = .2  # characteristic length [m]
res = 1000  # figure resolution
minval = 0.  # figure minimum value
maxval = 8.5  # figure maximum value
fac = 3  # animation duration multiplier

plt.figure(figsize=(10, 5))
lines = open(folder + 'cc.txt').readlines()
points = []
for line in lines:
    ccx = float(line.split()[0])
    ccy = float(line.split()[1])
    points.append([ccx, ccy])
grid_x, grid_y = np.mgrid[0:L:complex(0, res*L), 0:H:complex(0, res*H)]

nt = int(t/dt) + 1
for n in range(interval, nt, interval):
    lines = open(folder + variable + '{:010d}.txt'.format(n)).readlines()
    values = []
    for line in lines:  
        values.append(float(line.split()[0]))
    grid_z = griddata(points, values, (grid_x, grid_y), method='linear')

    plt.imshow(grid_z.T, extent=(0,L,0,H), origin='lower', vmin=minval, vmax=maxval, cmap='jet')
    plt.colorbar(orientation='horizontal', label='Speed [m/s]', shrink=0.5)
    plt.title('Re = %d, t = %.3f s' % (U*D/nu, n*dt), fontweight='bold')
    plt.xlabel('Channel Length [m]')
    plt.ylabel('Channel Height [m]')
    plt.axis('equal')
    plt.savefig(folder + variable + '{:010d}.png'.format(n))
    plt.clf()

    print('%d of %d (%.3f%%)' % (n, nt, n/nt*100))

image = []
for n in range(interval, nt, interval):
    image.append(folder + variable + '{:010d}.png'.format(n))
clip = mpy.ImageSequenceClip(image, fps=int(len(image)/t/fac))
clip.write_videofile(folder + 'animation.mp4')
