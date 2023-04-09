import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import moviepy.editor as mpy

mesh = 'cylinder.su2'  # mesh file
check = False  # check mesh

data = 'data/'  # data folder
interval = 50  # data interval

res = 1000  # figure resolution
minval = 0.  # figure minimum value
maxval = 8.5  # figure maximum value
fac = 3  # animation duration multiplier

U = 5.  # inlet x-velocity [m/s]
V = 0.  # inlet y-velocity [m/s]
nu = .01  # kinematic viscosity [m2/s]
t = 5.  # time [s]
dt = 1e-4  # time step [s]

L = 3.  # domain length [m]
H = 1.  # domain height [m]
D = .2  # characteristic length [m]

# program start

lines = open(mesh, 'r').readlines()

nodes = []
faces = []
cells = []

inlet_nodes = []
outlet_nodes = []
wall_nodes = []
body_nodes = []

cell_cells = []
cell_faces = []
face_cells = []
node_cells = []
n_node_cells = []

Vc = []  # cell volume
cc = []  # cell centre
signs = []  # cell face sign
Sf = []  # face area
nf = []  # face normal vector
tf = []  # face tangent vector
cf = []  # face centre
df = []  # cell-to-cell normal distance
ef = []  # cell-to-cell tangent distance
lf = []  # cell-to-cell distance
wf = []  # cell-to-face weighting factor
wn = []  # cell-to-node weighting factor

# cells

n_cells = int(lines[1].split()[1])
for i in range(n_cells):
    n1 = int(lines[2 + i].split()[1])
    n2 = int(lines[2 + i].split()[2])
    n3 = int(lines[2 + i].split()[3])
    cells.append([n1, n2, n3])

# nodes

n_nodes = int(lines[2 + n_cells].split()[1])
for i in range(n_nodes):
    nx = float(lines[3 + n_cells + i].split()[0])
    ny = float(lines[3 + n_cells + i].split()[1])
    nodes.append([nx, ny])

# faces

index = 5 + n_cells + n_nodes
n_inlet_faces = int(lines[index].split()[1])
for i in range(n_inlet_faces):
    n1 = int(lines[index + 1 + i].split()[1])
    n2 = int(lines[index + 1 + i].split()[2])
    faces.append([n1, n2])

    if n1 not in inlet_nodes:
        inlet_nodes.append(n1)
    if n2 not in inlet_nodes:
        inlet_nodes.append(n2)

index += 2 + n_inlet_faces
n_outlet_faces = int(lines[index].split()[1])
for i in range(n_outlet_faces):
    n1 = int(lines[index + 1 + i].split()[1])
    n2 = int(lines[index + 1 + i].split()[2])
    faces.append([n1, n2])

    if n1 not in outlet_nodes:
        outlet_nodes.append(n1)
    if n2 not in outlet_nodes:
        outlet_nodes.append(n2)

index += 2 + n_outlet_faces
n_wall_faces = int(lines[index].split()[1])
for i in range(n_wall_faces):
    n1 = int(lines[index + 1 + i].split()[1])
    n2 = int(lines[index + 1 + i].split()[2])
    faces.append([n1, n2])

    if n1 not in wall_nodes:
        wall_nodes.append(n1)
    if n2 not in wall_nodes:
        wall_nodes.append(n2)

index += 2 + n_wall_faces
n_body_faces = int(lines[index].split()[1])
for i in range(n_body_faces):
    n1 = int(lines[index + 1 + i].split()[2])  # swap
    n2 = int(lines[index + 1 + i].split()[1])  # swap
    faces.append([n1, n2])

    if n1 not in body_nodes:
        body_nodes.append(n1)
    if n2 not in body_nodes:
        body_nodes.append(n2)

for i in range(n_cells):
    n1 = cells[i][0]
    n2 = cells[i][1]
    n3 = cells[i][2]

    if [n1, n2] not in faces and [n2, n1] not in faces:
        faces.append([n1, n2])
    if [n2, n3] not in faces and [n3, n2] not in faces:
        faces.append([n2, n3])
    if [n3, n1] not in faces and [n1, n3] not in faces:
        faces.append([n3, n1])
n_faces = len(faces)

start = 0
end = n_inlet_faces
inlet_faces = list(range(start, end))

start += n_inlet_faces
end += n_outlet_faces
outlet_faces = list(range(start, end))

start += n_outlet_faces
end += n_wall_faces
wall_faces = list(range(start, end))

start += n_wall_faces
end += n_body_faces
body_faces = list(range(start, end))

start += n_body_faces
end = n_faces
interior_faces = list(range(start,end))
n_interior_faces = len(interior_faces)

boundary_faces = inlet_faces + outlet_faces + wall_faces + body_faces
n_boundary_faces = len(boundary_faces)

# cell_faces

for i in range(n_cells):
    n1 = cells[i][0]
    n2 = cells[i][1]
    n3 = cells[i][2]

    cell_faces_list = []
    for j in range(n_faces):
        if [n1, n2] == faces[j] or [n2, n1] == faces[j]:
            cell_faces_list.append(j)
        if [n2, n3] == faces[j] or [n3, n2] == faces[j]:
            cell_faces_list.append(j)
        if [n3, n1] == faces[j] or [n1, n3] == faces[j]:
            cell_faces_list.append(j)
    cell_faces.append(cell_faces_list)

# face_cells

for i in range(n_faces):
    face_cells_list = []
    for j in range(n_cells):
        if i in cell_faces[j]:
            face_cells_list.append(j)
    if len(face_cells_list) == 1:  # boundary face
        face_cells.append([face_cells_list[0], face_cells_list[0]])
    if len(face_cells_list) == 2:  # interior face
        face_cells.append([face_cells_list[0], face_cells_list[1]])

# Vc, cc, signs

for i in range(n_cells):
    n1 = cells[i][0]
    n2 = cells[i][1]
    n3 = cells[i][2]
    n1x = nodes[n1][0]
    n1y = nodes[n1][1]
    n2x = nodes[n2][0]
    n2y = nodes[n2][1]
    n3x = nodes[n3][0]
    n3y = nodes[n3][1]

    Vc.append(.5*abs((n1x*(n2y - n3y) + n2x*(n3y - n1y) + n3x*(n1y - n2y))))
    cc.append([(n1x + n2x + n3x)/3, (n1y + n2y + n3y)/3])

    signs_list = []
    cell_cells_list = []
    for j in range(3):
        fj = cell_faces[i][j]
        if i == face_cells[fj][0]:
            signs_list.append(1)
            cell_cells_list.append(face_cells[fj][1])
        else:
            signs_list.append(-1)
            cell_cells_list.append(face_cells[fj][0])
    signs.append(signs_list)
    cell_cells.append(cell_cells_list)

# Sf, nf, tf, cf, df, ef, lf, wf

for i in range(n_faces):
    n1 = faces[i][0]
    n2 = faces[i][1]
    n1x = nodes[n1][0]
    n1y = nodes[n1][1]
    n2x = nodes[n2][0]
    n2y = nodes[n2][1]

    Sf.append(m.sqrt((n2x - n1x)**2 + (n2y - n1y)**2))
    nf.append([(n2y - n1y)/Sf[i], -(n2x - n1x)/Sf[i]])
    tf.append([(n2x - n1x)/Sf[i], (n2y - n1y)/Sf[i]])
    cf.append([.5*(n1x + n2x), .5*(n1y + n2y)])

    c1 = face_cells[i][0]
    c2 = face_cells[i][1]
    cc1x = cc[c1][0]
    cc1y = cc[c1][1]
    cc2x = cc[c2][0]
    cc2y = cc[c2][1]
    cfx = cf[i][0]
    cfy = cf[i][1]

    if c1 == c2:  # boundary face
        lf.append([cfx - cc1x, cfy - cc1y])
        df.append(abs(np.dot(np.array(lf[i]), np.array(nf[i]))))
        ef.append(abs(np.dot(np.array(lf[i]), np.array(tf[i]))))
    else:  # interior face
        lf.append([cc2x - cc1x, cc2y - cc1y])
        df.append(abs(np.dot(np.array(lf[i]), np.array(nf[i]))))
        ef.append(abs(np.dot(np.array(lf[i]), np.array(tf[i]))))
    
    d1 = 1/m.sqrt((cc1x - cfx)**2 + (cc1y - cfy)**2)  # reciprocal
    d2 = 1/m.sqrt((cc2x - cfx)**2 + (cc2y - cfy)**2)  # reciprocal
    
    wf.append(d1/(d1 + d2))

# wn

for i in range(n_nodes):
    nx = nodes[i][0]
    ny = nodes[i][1]

    dj = []
    node_cells_list = []
    n_node_cells_value = 0
    for j in range(n_cells):
        if i in cells[j]:
            ccx = cc[j][0]
            ccy = cc[j][1]

            dj.append(1/m.sqrt((ccx - nx)**2 + (ccy - ny)**2))  # reciprocal

            node_cells_list.append(j)
            n_node_cells_value += 1
    node_cells.append(node_cells_list)
    n_node_cells.append(n_node_cells_value)

    wn_list = []
    for j in range(n_node_cells[i]):
        wn_list.append(dj[j]/sum(dj))
    wn.append(wn_list)

# check mesh

if check:
    ni = 0  # check node
    fi = 0  # check face
    ci = 0  # check cell

    plt.figure(figsize=(20, 10))

    for i in range(n_faces):
        n1 = faces[i][0]
        n2 = faces[i][1]
        n1x = nodes[n1][0]
        n1y = nodes[n1][1]
        n2x = nodes[n2][0]
        n2y = nodes[n2][1]

        cfx = cf[i][0]
        cfy = cf[i][1]

        plt.text(cfx, cfy, i, horizontalalignment='center', verticalalignment='center')

        if i in boundary_faces:
            plt.plot([n1x, n2x], [n1y, n2y], c='r', linewidth=2.)
        else:
            plt.plot([n1x, n2x], [n1y, n2y], c='k', linewidth=.5)

        if i == fi:
            c1 = face_cells[i][0]
            c2 = face_cells[i][1]
            cc1x = cc[c1][0]
            cc1y = cc[c1][1]
            cc2x = cc[c2][0]
            cc2y = cc[c2][1]

            plt.text(n1x, n1y, 1, c='b', horizontalalignment='center', verticalalignment='center')
            plt.text(n2x, n2y, 2, c='b', horizontalalignment='center', verticalalignment='center')
            plt.text(cc1x, cc1y, 1, c='r', horizontalalignment='center', verticalalignment='center')
            plt.text(cc2x, cc2y, 2, c='r', horizontalalignment='center', verticalalignment='center')
            plt.arrow(cfx, cfy, .1*nf[i][0], .1*nf[i][1], color='g', head_width=.02)
            plt.arrow(cfx, cfy, .1*tf[i][0], .1*tf[i][1], color='y', head_width=.02)
            
            plt.arrow(cc1x, cc1y, df[i]*nf[i][0], df[i]*nf[i][1], color='c', head_width=.02)
            plt.arrow(cc1x, cc1y, ef[i]*tf[i][0], ef[i]*tf[i][1], color='m', head_width=.02)
            plt.arrow(cc1x, cc1y, lf[i][0], lf[i][1], color='k', head_width=.02)

            print('Face', fi, c1, c2, wf[i])

    for i in range(n_cells):
        ccx = cc[i][0]
        ccy = cc[i][1]

        plt.text(ccx, ccy, i, horizontalalignment='center', verticalalignment='center')

        if i == ci:
            for j in range(3):
                fj = cell_faces[ci][j]
                cj = cell_cells[ci][j]
                c1 = face_cells[fj][0]
                c2 = face_cells[fj][1]

                print('Cell', ci, fj, cj, c1, c2, signs[ci][j])
    
    for i in range(n_nodes):
        nx = nodes[i][0]
        ny = nodes[i][1]

        plt.text(nx, ny, i, horizontalalignment='center', verticalalignment='center')

        if i == ni:
            for j in range(n_node_cells[i]):
                cj = node_cells[i][j]
                ccjx = cc[cj][0]
                ccjy = cc[cj][1]

                plt.scatter(ccjx, ccjy)
                print('Node', ni, cj, wn[i][j], n_node_cells[i])

    plt.axis('equal')
    plt.grid('on')
    plt.show()

u = np.zeros(n_cells)
us = np.zeros(n_cells)
un = np.zeros(n_nodes)

v = np.zeros(n_cells)
vs = np.zeros(n_cells)
vn = np.zeros(n_nodes)

p = np.zeros(n_cells)
pf = np.zeros(n_faces)

Jcu = np.zeros(n_faces)
Jcv = np.zeros(n_faces)
Jdu = np.zeros(n_faces)
Jdv = np.zeros(n_faces)

aC = np.zeros(n_cells)
aF = np.zeros((n_cells, 3))
bC = np.zeros(n_cells)

# pressure coefficients

for i in range(n_cells):
    aC[i] = 0.
    for j in range(3):
        fj = cell_faces[i][j]

        if fj in inlet_faces:
            aF[i, j] = 0.
            aC[i] += 0.
        elif fj in outlet_faces:
            aF[i, j] = 0.
            aC[i] += -Sf[fj]/df[fj]
        elif fj in wall_faces or fj in body_faces:
            aF[i, j] = 0.
            aC[i] += 0.
        else:  # interior face
            aF[i, j] = Sf[fj]/df[fj]
            aC[i] += -Sf[fj]/df[fj]

# time loop

Nt = int(t/dt)
for nt in range(Nt):
    
    # projection step

    for i in range(n_nodes):
        un[i] = 0.
        vn[i] = 0.
        for j in range(n_node_cells[i]):
            cj = node_cells[i][j]
            
            un[i] += wn[i][j]*u[cj]
            vn[i] += wn[i][j]*v[cj]
    for i in inlet_nodes:
        un[i] = U
        vn[i] = V
    for i in wall_nodes:
        un[i] = 0.
        vn[i] = 0.
    for i in body_nodes:
        un[i] = 0.
        vn[i] = 0.

    for i in range(n_interior_faces):
        fi = interior_faces[i]
        c1 = face_cells[fi][0]
        c2 = face_cells[fi][1]
        n1 = faces[fi][0]
        n2 = faces[fi][1]

        Sfi = Sf[fi]
        dfi = df[fi]
        efi = ef[fi]

        Jdu[fi] = nu*((u[c2] - u[c1])*Sfi - (un[n2] - un[n1])*efi)/dfi
        Jdv[fi] = nu*((v[c2] - v[c1])*Sfi - (vn[n2] - vn[n1])*efi)/dfi

        wfi = wf[fi]
        ufi = wfi*u[c1] + (1 - wfi)*u[c2]
        vfi = wfi*v[c1] + (1 - wfi)*v[c2]
        Vfi = np.array([ufi, vfi])
        nfi = np.array(nf[fi])
        mfi = np.dot(Vfi, nfi)

        Jcu[fi] = mfi*Sfi*ufi
        Jcv[fi] = mfi*Sfi*vfi

    for i in range(n_boundary_faces):
        fi = boundary_faces[i]
        ci = face_cells[fi][0]
        
        if fi in inlet_faces:
            Sfi = Sf[fi]
            dfi = df[fi]

            Jdu[fi] = nu*(U - u[ci])*Sfi/dfi
            Jdv[fi] = nu*(V - v[ci])*Sfi/dfi

            Vfi = np.array([U, V])
            nfi = np.array(nf[fi])
            mfi = np.dot(Vfi, nfi)

            Jcu[fi] = mfi*Sfi*U
            Jcv[fi] = mfi*Sfi*V
        
        if fi in outlet_faces:
            Sfi = Sf[fi]

            Jdu[fi] = 0.
            Jdv[fi] = 0.

            Vfi = np.array([u[ci], v[ci]])
            nfi = np.array(nf[fi])
            mfi = np.dot(Vfi, nfi)

            Jcu[fi] = mfi*Sfi*u[ci]
            Jcv[fi] = mfi*Sfi*v[ci]
        
        if fi in wall_faces or fi in body_faces:
            Sfi = Sf[fi]
            dfi = df[fi]

            Jdu[fi] = -nu*u[ci]*Sfi/dfi
            Jdv[fi] = -nu*v[ci]*Sfi/dfi

            Jcu[fi] = 0.
            Jcv[fi] = 0.

    for i in range(n_cells):
        us[i] = u[i]
        vs[i] = v[i]
        for j in range(3):
            fj = cell_faces[i][j]

            sign = signs[i][j]

            us[i] += sign*dt/Vc[i]*(Jdu[fj] - Jcu[fj])
            vs[i] += sign*dt/Vc[i]*(Jdv[fj] - Jcv[fj])
    
    # correction step

    for i in range(n_cells):
        bC[i] = 0.
        for j in range(3):
            fj = cell_faces[i][j]
            c1 = face_cells[fj][0]
            c2 = face_cells[fj][1]

            if fj in wall_faces or fj in body_faces:
                bC[i] += 0.
            else:
                Sfj = Sf[fj]
                sign = signs[i][j]

                wfj = wf[fj]
                usfj = wfj*us[c1] + (1 - wfj)*us[c2]
                vsfj = wfj*vs[c1] + (1 - wfj)*vs[c2]
                Vsfj = np.array([usfj, vsfj])
                nfj = np.array(nf[fj])
                msfj = np.dot(Vsfj, nfj)

                bC[i] += sign*msfj*Sfj/dt

    for k in range(5):  # Gauss-Seidel
        for i in range(n_cells):
            value = 0.
            for j in range(3):
                cj = cell_cells[i][j]

                value += aF[i, j]*p[cj]
            p[i] = (bC[i] - value)/aC[i]

    for i in range(n_faces):
        c1 = face_cells[i][0]
        c2 = face_cells[i][1]
        
        pf[i] = wf[i]*p[c1] + (1 - wf[i])*p[c2]
    for i in outlet_faces:
        pf[i] = 0.

    # update step

    for i in range(n_cells):
        u[i] = us[i]
        v[i] = vs[i]
        for j in range(3):
            fj = cell_faces[i][j]

            Sfj = Sf[fj]
            sign = signs[i][j]
            nfj = np.array(nf[fj])

            u[i] += -dt/Vc[i]*pf[fj]*sign*nfj[0]*Sfj
            v[i] += -dt/Vc[i]*pf[fj]*sign*nfj[1]*Sfj

    print('%d of %d (%.3f%%)' % (nt, Nt, nt/Nt*100))
    print(np.max(u), np.min(u))

    if nt % interval == 0:
        np.save(data + 'u{:d}'.format(nt), u)
        np.save(data + 'v{:d}'.format(nt), v)

plt.figure(figsize=(10, 5))
points = np.zeros((n_cells, 2))
for i in range(n_cells):
    points[i, 0] = cc[i][0]
    points[i, 1] = cc[i][1]
grid_x, grid_y = np.mgrid[0:L:complex(0, res*L), 0:H:complex(0, res*H)]

for nt in range(0, Nt, interval):
    u = np.load(data + 'u{:d}.npy'.format(nt))
    v = np.load(data + 'v{:d}.npy'.format(nt))
    values = np.sqrt(u*u + v*v)
    grid_z = griddata(points, values, (grid_x, grid_y), method='linear')
    
    plt.imshow(grid_z.T, extent=(0,L,0,H), origin='lower', vmin=minval, vmax=maxval, cmap='jet')
    plt.colorbar(orientation='horizontal', label='Speed [m/s]', shrink=0.5)
    plt.title('Re = %d, t = %.3f s' % (U*D/nu, nt*dt), fontweight='bold')
    plt.xlabel('Channel Length [m]')
    plt.ylabel('Channel Height [m]')
    plt.axis('equal')
    plt.savefig(data + 'figure{:d}.png'.format(nt))
    plt.clf()

    print('%d of %d (%.3f%%)' % (nt, Nt, nt/Nt*100))

image = []
for nt in range(0, Nt, interval):
    image.append(data + 'figure{:d}.png'.format(nt))
clip = mpy.ImageSequenceClip(image, fps=int(len(image)/t/fac))
clip.write_videofile(data + 'animation.mp4')
