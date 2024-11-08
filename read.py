import numpy as np
import matplotlib.pyplot as plt
import moviepy.editor as mpy
from scipy.interpolate import griddata

xb = [1., 4.]           # grid x-bounds
yb = [1.5, 3.5]         # grid y-bounds
vb = [0., 3.5]          # velocity bounds
wb = [-150., 150.]      # vorticity bounds
fs = [8, 8]             # figure size
res = 200               # grid resolution
fps = 100               # frames per second
full = True             # full plots

# read configuration

with open('config.txt', 'r') as file:
    for line in file:
        if line.startswith('Re='):
            Re = float(line.split()[1])
        elif line.startswith('nt='):
            nt = int(line.split()[1])
        elif line.startswith('nw='):
            nw = int(line.split()[1])

# field animation

print('Generating animation...')
nxy = np.loadtxt('data/nxy.txt')
cxy = np.loadtxt('data/cxy.txt')
x, y = cxy[:, 0], cxy[:, 1]
xx, yy = np.mgrid[xb[0]:xb[1]:res*1j, yb[0]:yb[1]:res*1j]
dx = xx[1, 0] - xx[0, 0]
dy = yy[0, 1] - yy[0, 0]

within_box = (x > xb[0]) & (x < xb[1]) & (y > yb[0]) & (y < yb[1])
indices = np.where(within_box)[0]

frames = []
w = np.zeros((res, res)) # vorticity
for i in range(nw, nt + nw, nw):
    data = np.loadtxt(f'data/Wc_{i:010d}.txt')
    
    # compute velocity
    uxy = np.sqrt(data[indices, 0]**2 + data[indices, 1]**2)
    v = griddata((x[indices], y[indices]), uxy, (xx, yy), method='linear')
    # compute vorticity
    ux, uy = data[indices, 0], data[indices, 1]
    uxg = griddata((x[indices], y[indices]), ux, (xx, yy), method='linear')
    uyg = griddata((x[indices], y[indices]), uy, (xx, yy), method='linear')
    dvdx = (uyg[2:, :] - uyg[:-2, :])/dx
    dudy = (uxg[:, 2:] - uxg[:, :-2])/dy
    w[1:-1, 1:-1] = dvdx[:, 1:-1] - dudy[1:-1, :]

    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=fs)
    # plot velocity
    cf1 = ax1.contourf(xx, yy, v, \
    levels=np.linspace(vb[0], vb[1], res), cmap='jet', extend='both')
    ax1.fill(nxy[:, 0], nxy[:, 1])
    ax1.set_aspect('equal')
    # plot vorticity
    cf2 = ax2.contourf(xx, yy, w, \
    levels=np.linspace(wb[0], wb[1], res), cmap='jet', extend='both')
    ax2.fill(nxy[:, 0], nxy[:, 1])
    ax2.set_aspect('equal')

    if full:
        text1 = 'Speed (Re = {:.0f})'.format(Re)
        ax1.text(.025, .025, text1, fontweight='bold', \
        transform=ax1.transAxes, va='bottom', ha='left')
        cb1 = fig.colorbar(cf1, ax=ax1)
        cb1.set_ticks(np.linspace(vb[0], vb[1], 11))

        text2 = 'Vorticity (Re = {:.0f})'.format(Re)
        ax2.text(0.025, .025, text2, fontweight='bold', \
        transform=ax2.transAxes, va='bottom', ha='left')
        cb2 = fig.colorbar(cf2, ax=ax2)
        cb2.set_ticks(np.linspace(wb[0], wb[1], 11))

        fig.tight_layout()
        fig.savefig(f'data/frame_{i:010d}.png')
        fig.clf()
    else:
        ax1.axis('off')
        ax2.axis('off')

        fig.tight_layout()
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.savefig(f'data/frame_{i:010d}.png', bbox_inches='tight', pad_inches=0)
        fig.clf()

    frames.append(f'data/frame_{i:010d}.png')
    print('Saving frame %d of %d (%.3f%%) ... Done!' % (i, nt, i/nt*100.))

animation = mpy.ImageSequenceClip(frames, fps=fps)
animation.write_videofile('data/animation.mp4')
print('Done!')
