import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import griddata

xb = [1., 4.]       # grid x-bounds
yb = [1., 2.]       # grid y-bounds
fsz = [12, 4]       # figure size
res = 200           # grid resolution
fps = 100           # frames per second

# read configuration

with open("config.txt", "r") as file:
    for line in file:
        if line.startswith("nt="):
            nt = int(line.split()[1])
        elif line.startswith("nw="):
            nw = int(line.split()[1])

# field animation

print('Generating animation...')
nxy = np.loadtxt('data/nxy.txt')
cxy = np.loadtxt('data/cxy.txt')
x, y = cxy[:, 0], cxy[:, 1]
within_box = (x > xb[0]) & (x < xb[1]) & (y > yb[0]) & (y < yb[1])
indices = np.where(within_box)[0]

grid_x, grid_y = np.mgrid[xb[0]:xb[1]:res*1j, yb[0]:yb[1]:res*1j]
data = np.loadtxt(f'data/Wc_{nw:010d}.txt')
values = np.sqrt(data[indices, 0]**2 + data[indices, 1]**2)
grid_values = griddata((x[indices], y[indices]), values, (grid_x, grid_y), method='linear')

fig, ax = plt.subplots(figsize=fsz)
ax.set_aspect('equal')
ax.set_title('Speed [m/s]', fontweight='bold')
contour = ax.contourf(grid_x, grid_y, grid_values, levels=res, cmap='jet')
cbar = plt.colorbar(contour, ax=ax)
ax.fill(nxy[:, 0], nxy[:, 1])

def update(frame):
    data = np.loadtxt(f'data/Wc_{frame:010d}.txt')
    values = np.sqrt(data[indices, 0]**2 + data[indices, 1]**2)
    grid_values = griddata((x[indices], y[indices]), values, (grid_x, grid_y), method='linear')
    
    ax.cla()
    ax.set_aspect('equal')
    ax.set_title('Speed [m/s]', fontweight='bold')
    contour = ax.contourf(grid_x, grid_y, grid_values, levels=res, cmap='jet')
    cbar.update_normal(contour)
    cbar.update_ticks()
    ax.fill(nxy[:, 0], nxy[:, 1])

    return contour

ani = FuncAnimation(fig, update, frames=range(nw, nt + nw, nw), blit=False)
ani.save('data/speed.mp4', writer='ffmpeg', fps=fps)
print('Done!')
plt.show()
