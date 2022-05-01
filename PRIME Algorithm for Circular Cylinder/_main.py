import math as m
import numpy as np
import random as rd
import matplotlib.pyplot as plt
import matplotlib.colors as clr

U = 1.0  # inlet velocity [m/s]
nu = 0.01  # kinematic viscosity [m2/s]
L = 0.2  # characteristic length [m]

itr = 5000  # number of iterations
prf = 1.0  # pressure relaxation factor

mesh = 'cylinder.su2'  # mesh file [.su2]
check_mesh = 0  # check mesh

###################
## program start ##
###################

print('Setting up mesh |..........| (0%)')

file = open(mesh, 'r')  # open file
lines = file.readlines()  # read lines

nodes = []  # nodes [[nx, ny], ...]
faces = []  # faces [[n1, n2], ...]
cells = []  # cells [[n1, n2, n3], ...]
n_nodes = 0  # number of nodes
n_faces = 0  # number of faces
n_cells = 0  # number of cells

wall_faces = []  # wall faces [f1, f2, ...]
object_faces = []  # object faces [f1, f2, ...]
inlet_faces = []  # inlet faces [f1, f2, ...]
outlet_faces = []  # outlet faces [f1, f2, ...]
boundary_faces = []  # boundary faces [f1, f2, ...]
interior_faces = []  # interior faces [f1, f2, ...]
n_wall_faces = 0  # number of wall faces
n_object_faces = 0  # number of object faces
n_inlet_faces = 0  # number of inlet faces
n_outlet_faces = 0  # number of outlet faces
n_boundary_faces = 0  # number of boundary faces
n_interior_faces = 0  # number of interior faces

boundary_cells = []  # boundary cells [c1, c2, ...]
interior_cells = []  # interior cells [c1, c2, ...]
n_boundary_cells = 0  # number of boundary cells
n_interior_cells = 0  # number of interior cells

cell_faces = []  # faces surrounding cell [[f1, f2, f3], ...]
cell_cells = []  # cells surrounding cell [[c1, c2, c3], ...]
face_cells = []  # cells surrounding face [[c1, c2], ...]

cc = []  # cell centres [[ccx, ccy], ...]
Vc = []  # cell volumes [Vc1, Vc2, ...]
nt = []  # cell faces normal/tangent signs [[nt1, nt2, nt3], ...]

fc = []  # face centres [[fcx, fcy], ...]
Sf = []  # face surface areas [Af1, Af2, ...]
nf = []  # face surface normals [[nfx, nfy], ...]
tf = []  # face surface tangents [[tfx, tfy], ...]
df = []  # face cells distances [df1, df2, ...]
ef = []  # face diffusion unit vectors [[efx, efy], ...]
Ef = []  # face diffusion vector lengths [Ef1, Ef2, ...]
Tf = []  # face cross diffusion vectors [[Tfx, Tfy], ...]
gf = []  # face weighting factors [gf1, gf2, ...]
hf = []  # face cells perpendicular distances [hf1, hf2, ...]

# nodes, faces, cells
for i in range(len(lines)):
    if lines[i].split()[0] == 'NPOIN=':
        n_nodes = int(lines[i].split()[1])
        for j in range(n_nodes):
            nx = float(lines[i + j + 1].split()[0])  # node x-coordinate
            ny = float(lines[i + j + 1].split()[1])  # node y-coordinate
            nodes.append([nx, ny])
        break

for i in range(len(lines)):
    if lines[i].split()[1].upper() == 'WALL':
        n_wall_faces = int(lines[i + 1].split()[1])
        for j in range(n_wall_faces):
            n1 = int(lines[i + j + 2].split()[1])  # face node 1
            n2 = int(lines[i + j + 2].split()[2])  # face node 2
            faces.append([n1, n2])
        break

for i in range(len(lines)):
    if lines[i].split()[1].upper() == 'OBJECT':
        n_object_faces = int(lines[i + 1].split()[1])
        for j in range(n_object_faces):
            n1 = int(lines[i + j + 2].split()[2])  # face node 1 (swap)
            n2 = int(lines[i + j + 2].split()[1])  # face node 2 (swap)
            faces.append([n1, n2])
        break

for i in range(len(lines)):
    if lines[i].split()[1].upper() == 'INLET':
        n_inlet_faces = int(lines[i + 1].split()[1])
        for j in range(n_inlet_faces):
            n1 = int(lines[i + j + 2].split()[1])  # face node 1
            n2 = int(lines[i + j + 2].split()[2])  # face node 2
            faces.append([n1, n2])
        break

for i in range(len(lines)):
    if lines[i].split()[1].upper() == 'OUTLET':
        n_outlet_faces = int(lines[i + 1].split()[1])
        for j in range(n_outlet_faces):
            n1 = int(lines[i + j + 2].split()[1])  # face node 1
            n2 = int(lines[i + j + 2].split()[2])  # face node 2
            faces.append([n1, n2])
        break

for i in range(len(lines)):
    if lines[i].split()[0] == 'NELEM=':
        n_cells = int(lines[i].split()[1])
        for j in range(n_cells):
            n1 = int(lines[i + j + 1].split()[1])  # cell node 1
            n2 = int(lines[i + j + 1].split()[2])  # cell node 2
            n3 = int(lines[i + j + 1].split()[3])  # cell node 3
            cells.append([n1, n2, n3])

            if [n1, n2] not in faces and [n2, n1] not in faces:
                faces.append([n1, n2])
            if [n2, n3] not in faces and [n3, n2] not in faces:
                faces.append([n2, n3])
            if [n3, n1] not in faces and [n1, n3] not in faces:
                faces.append([n3, n1])
        n_faces = len(faces)
        break

start = 0
end = n_wall_faces
wall_faces = list(range(start, end))

start += n_wall_faces
end += n_object_faces
object_faces = list(range(start, end))

start += n_object_faces
end += n_inlet_faces
inlet_faces = list(range(start, end))

start += n_inlet_faces
end += n_outlet_faces
outlet_faces = list(range(start, end))

start += n_outlet_faces
end = n_faces
interior_faces = list(range(start, end))
n_interior_faces = len(interior_faces)

boundary_faces = wall_faces + object_faces + inlet_faces + outlet_faces
n_boundary_faces = len(boundary_faces)

print('Setting up mesh |##........| (20%)')

# cell_faces
for i in range(n_cells):
    n1 = cells[i][0]  # cell node 1
    n2 = cells[i][1]  # cell node 2
    n3 = cells[i][2]  # cell node 3

    cell_faces_list = []
    for j in range(n_faces):
        if [n1, n2] == faces[j] or [n2, n1] == faces[j]:
            cell_faces_list.append(j)
        if [n2, n3] == faces[j] or [n3, n2] == faces[j]:
            cell_faces_list.append(j)
        if [n3, n1] == faces[j] or [n1, n3] == faces[j]:
            cell_faces_list.append(j)
    cell_faces.append(cell_faces_list)

print('Setting up mesh |####......| (40%)')

# face_cells
for i in range(n_faces):
    face_cells_list = []
    for j in range(n_cells):
        if i in cell_faces[j]:
            face_cells_list.append(j)
    if len(face_cells_list) == 1:  # boundary face
        face_cells.append([face_cells_list[0], face_cells_list[0]])
    if len(face_cells_list) == 2:  # interior face
        face_cells.append([face_cells_list[1], face_cells_list[0]])

print('Setting up mesh |######....| (60%)')

# cc, Vc, nt
for i in range(n_cells):
    n1 = cells[i][0]  # cell node 1
    n2 = cells[i][1]  # cell node 1
    n3 = cells[i][2]  # cell node 1

    n1x = nodes[n1][0]  # cell node 1 x-coordinate
    n1y = nodes[n1][1]  # cell node 1 y-coordinate
    n2x = nodes[n2][0]  # cell node 2 x-coordinate
    n2y = nodes[n2][1]  # cell node 2 y-coordinate
    n3x = nodes[n3][0]  # cell node 3 x-coordinate
    n3y = nodes[n3][1]  # cell node 3 y-coordinate

    ccx = (1/3)*(n1x + n2x + n3x)  # cell centre x-coordinate
    ccy = (1/3)*(n1y + n2y + n3y)  # cell centre y-coordinate
    cc.append([ccx, ccy])

    Vc.append((1/2)*abs((n1x*(n2y - n3y) + n2x*(n3y - n1y) + n3x*(n1y - n2y))))

    f1 = cell_faces[i][0]  # cell face 1
    f2 = cell_faces[i][1]  # cell face 2
    f3 = cell_faces[i][2]  # cell face 3

    nt_list = []
    cell_cells_list = []
    if i == face_cells[f1][0]:
        nt_list.append(1)
        cell_cells_list.append(face_cells[f1][1])
    else:
        nt_list.append(-1)
        cell_cells_list.append(face_cells[f1][0])

    if i == face_cells[f2][0]:
        nt_list.append(1)
        cell_cells_list.append(face_cells[f2][1])
    else:
        nt_list.append(-1)
        cell_cells_list.append(face_cells[f2][0])

    if i == face_cells[f3][0]:
        nt_list.append(1)
        cell_cells_list.append(face_cells[f3][1])
    else:
        nt_list.append(-1)
        cell_cells_list.append(face_cells[f3][0])
    nt.append(nt_list)
    cell_cells.append(cell_cells_list)

    if f1 in boundary_faces or f2 in boundary_faces or f3 in boundary_faces:
        boundary_cells.append(i)
    else:
        interior_cells.append(i)
n_boundary_cells = len(boundary_cells)
n_interior_cells = len(interior_cells)

print('Setting up mesh |########..| (80%)')

# fc, Sf, nf, tf, df, ef, Ef, Tf, gf, hf
for i in range(n_faces):
    n1 = faces[i][0]  # face node 1
    n2 = faces[i][1]  # face node 2

    n1x = nodes[n1][0]  # face node 1 x-coordinate
    n1y = nodes[n1][1]  # face node 1 y-coordinate
    n2x = nodes[n2][0]  # face node 2 x-coordinate
    n2y = nodes[n2][1]  # face node 2 y-coordinate

    fcx = (1/2)*(n1x + n2x)  # face centre x-coordinate
    fcy = (1/2)*(n1y + n2y)  # face centre y-coordinate
    fc.append([fcx, fcy])

    Sf.append(m.sqrt((n2x - n1x)**2 + (n2y - n1y)**2))
    nf.append([(n2y - n1y)/Sf[i], -(n2x - n1x)/Sf[i]])
    tf.append([(n2x - n1x)/Sf[i], (n2y - n1y)/Sf[i]])

    c1 = face_cells[i][0]  # face cell 1
    c2 = face_cells[i][1]  # face cell 2

    cc1x = cc[c1][0]  # cell 1 centre x-coordinate
    cc1y = cc[c1][1]  # cell 1 centre y-coordinate
    cc2x = cc[c2][0]  # cell 2 centre x-coordinate
    cc2y = cc[c2][1]  # cell 2 centre y-coordinate

    if c1 == c2:  # boundary faces
        df.append(m.sqrt((fcx - cc1x)**2 + (fcy - cc1y)**2))
        ef.append([(fcx - cc1x)/df[i], (fcy - cc1y)/df[i]])
        Ef.append(np.dot(Sf[i]*np.array(nf[i]), Sf[i]*np.array(nf[i]))/np.dot(Sf[i]*np.array(nf[i]), np.array(ef[i])))
        Tf.append([Sf[i]*nf[i][0] - Ef[i]*ef[i][0], Sf[i]*nf[i][1] - Ef[i]*ef[i][1]])
        hf.append(np.dot(df[i]*np.array(ef[i]), np.array(nf[i])))
    else:  # interior faces
        df.append(m.sqrt((cc2x - cc1x)**2 + (cc2y - cc1y)**2))
        ef.append([(cc2x - cc1x)/df[i], (cc2y - cc1y)/df[i]])
        Ef.append(np.dot(Sf[i]*np.array(nf[i]), Sf[i]*np.array(nf[i]))/np.dot(Sf[i]*np.array(nf[i]), np.array(ef[i])))
        Tf.append([Sf[i]*nf[i][0] - Ef[i]*ef[i][0], Sf[i]*nf[i][1] - Ef[i]*ef[i][1]])
        hf.append(np.dot(df[i]*np.array(ef[i]), np.array(nf[i])))

    Vc1 = Vc[c1]  # cell C
    Vc2 = Vc[c2]  # cell F
    gf.append(Vc1/(Vc1 + Vc2))

print('Setting up mesh |##########| (100%)')

if check_mesh:
    fi = rd.randint(0, n_faces)  # inspected face
    ci = rd.randint(0, n_cells)  # inspected cell

    for i in range(n_faces):
        n1x = nodes[faces[i][0]][0]  # face node 1 x-coordinate
        n1y = nodes[faces[i][0]][1]  # face node 1 y-coordinate
        n2x = nodes[faces[i][1]][0]  # face node 2 x-coordinate
        n2y = nodes[faces[i][1]][1]  # face node 2 y-coordinate
        plt.plot([n1x, n2x], [n1y, n2y], c='grey', lw=0.5)

    for i in range(n_cells):
        ccx = cc[i][0]  # cell centre x-coordinate
        ccy = cc[i][1]  # cell centre y-coordinate
        plt.text(ccx, ccy, i, c='b', horizontalalignment='center', verticalalignment='center')

    for i in range(n_faces):
        fcx = fc[i][0]  # face centre x-coordinate
        fcy = fc[i][1]  # face centre y-coordinate
        plt.text(fcx, fcy, i, c='r', horizontalalignment='center', verticalalignment='center')

    for i in range(n_nodes):
        nx = nodes[i][0]  # node x-coordinate
        ny = nodes[i][1]  # node y-coordinate
        plt.text(nx, ny, i, c='k', horizontalalignment='center', verticalalignment='center')

    # for i in interior_cells:
    #     ccx = cc[i][0]  # cell centre x-coordinate
    #     ccy = cc[i][1]  # cell centre y-coordinate
    #     plt.scatter(ccx, ccy, c='b', zorder=10)

    # for i in interior_faces:
    #     n1x = nodes[faces[i][0]][0]  # face node 1 x-coordinate
    #     n1y = nodes[faces[i][0]][1]  # face node 1 y-coordinate
    #     n2x = nodes[faces[i][1]][0]  # face node 2 x-coordinate
    #     n2y = nodes[faces[i][1]][1]  # face node 2 y-coordinate
    #     plt.plot([n1x, n2x], [n1y, n2y], c='r', lw=1.0)

    print('Face', fi)
    print('Nodes', faces[fi])
    print('Cells', face_cells[fi])
    plt.arrow(fc[fi][0], fc[fi][1], 0.1*nf[fi][0], 0.1*nf[fi][1], color='green', head_width=0.02)
    plt.arrow(fc[fi][0], fc[fi][1], 0.1*tf[fi][0], 0.1*tf[fi][1], color='orange', head_width=0.02)

    print('Cell', ci)
    print('Faces', cell_faces[ci])
    print('Cells', cell_cells[ci])
    print('Signs', nt[ci])
    for j in range(3):
        fj = cell_faces[ci][j]
        sign = nt[ci][j]
        plt.arrow(fc[fj][0], fc[fj][1], sign*Sf[fj]*nf[fj][0], sign*Sf[fj]*nf[fj][1], color='green', head_width=0.02)
        plt.arrow(fc[fj][0], fc[fj][1], sign*Sf[fj]*tf[fj][0], sign*Sf[fj]*tf[fj][1], color='orange', head_width=0.02)

        plt.arrow(cc[ci][0], cc[ci][1], sign*df[fj]*ef[fj][0], sign*df[fj]*ef[fj][1], color='grey', head_width=0.02)
        plt.arrow(cc[ci][0], cc[ci][1], sign*hf[fj]*nf[fj][0], sign*hf[fj]*nf[fj][1], color='black', head_width=0.02)
        plt.arrow(fc[fj][0], fc[fj][1], sign*Ef[fj]*ef[fj][0], sign*Ef[fj]*ef[fj][1], color='cyan', head_width=0.02)
        plt.arrow(fc[fj][0] + sign*Ef[fj]*ef[fj][0], fc[fj][1] + sign*Ef[fj]*ef[fj][1], sign*Tf[fj][0], sign*Tf[fj][1], color='magenta', head_width=0.02)

    plt.axis('equal')
    plt.show()

u = np.zeros(n_cells)  # cell u-velocity
us = np.zeros(n_cells)  # cell intermediate u-velocity
uf = np.zeros(n_faces)  # face u-velocity

v = np.zeros(n_cells)  # cell v-velocity
vs = np.zeros(n_cells)  # cell intermediate v-velocity
vf = np.zeros(n_faces)  # face v-velocity

p = np.zeros(n_cells)  # cell pressure
pf = np.zeros(n_faces)  # face pressure
ppf = np.zeros(n_faces)  # face pressure correction

mf = np.zeros(n_faces)  # face mass flow rate
msf = np.zeros(n_faces)  # face intermediate mass flow rate

grad_u = np.zeros((n_cells, 2))  # cell u-velocity gradient
grad_uf = np.zeros((n_faces, 2))  # face u-velocity gradient
grad_v = np.zeros((n_cells, 2))  # cell v-velocity gradient
grad_vf = np.zeros((n_faces, 2))  # face v-velocity gradient
grad_p = np.zeros((n_cells, 2))  # cell pressure gradient
grad_pf = np.zeros((n_faces, 2))  # face pressure gradient
grad_pp = np.zeros((n_cells, 2))  # cell pressure correction gradient

aCu = np.zeros(n_cells)  # current cell u-velocity coefficient
aFu = np.zeros((n_cells, 3))  # neighbour cell u-velocity coefficient
bCu = np.zeros(n_cells)  # current cell u-velocity source term
aCv = np.zeros(n_cells)  # current cell v-velocity coefficient
aFv = np.zeros((n_cells, 3))  # neighbour cell v-velocity coefficient
bCv = np.zeros(n_cells)  # current cell v-velocity source term

DCu = np.zeros(n_cells)
DCv = np.zeros(n_cells)
Dfu = np.zeros(n_faces)
Dfv = np.zeros(n_faces)
Epf = np.zeros(n_faces)

aCp = np.zeros(n_cells)  # current cell pressure correction coefficient
aFp = np.zeros((n_cells, 3))  # neighbour cell pressure correction coefficient
bCp = np.zeros(n_cells)  # current cell pressure correction source term

mf_in = 0.0
for fi in inlet_faces:
    mf[fi] = -Sf[fi]*U
    msf[fi] = -Sf[fi]*U
    mf_in += Sf[fi]*U

err = np.zeros(itr)
for n in range(itr):

    #######################
    ## velocity gradient ##
    #######################

    # face velocity
    for fi in inlet_faces:
        uf[fi] = U
        # vf[fi] = V
    for fi in interior_faces + outlet_faces:
        c1 = face_cells[fi][0]  # cell C
        c2 = face_cells[fi][1]  # cell F
        gfi = gf[fi]
        uf[fi] = gfi*u[c2] + (1 - gfi)*u[c1]
        vf[fi] = gfi*v[c2] + (1 - gfi)*v[c1]
    # cell velocity gradient
    for ci in range(n_cells):
        Vci = Vc[ci]
        grad_u[ci] = 0.0
        grad_v[ci] = 0.0
        for j in range(3):
            fj = cell_faces[ci][j]
            Sfj = Sf[fj]
            nfj = nt[ci][j]*np.array(nf[fj])
            grad_u[ci] += (1/Vci)*uf[fj]*Sfj*nfj
            grad_v[ci] += (1/Vci)*vf[fj]*Sfj*nfj
    # face velocity gradient
    for fi in interior_faces:
        c1 = face_cells[fi][0]  # cell C
        c2 = face_cells[fi][1]  # cell F
        gfi = gf[fi]
        grad_uf_bar = gfi*grad_u[c2] + (1 - gfi)*grad_u[c1]
        grad_vf_bar = gfi*grad_v[c2] + (1 - gfi)*grad_v[c1]

        dfi = df[fi]
        efi = np.array(ef[fi])
        grad_uf[fi] = grad_uf_bar + ((u[c2] - u[c1])/dfi - np.dot(grad_uf_bar, efi))*efi
        grad_vf[fi] = grad_vf_bar + ((v[c2] - v[c1])/dfi - np.dot(grad_vf_bar, efi))*efi
    for fi in inlet_faces + wall_faces + object_faces:
        c1 = face_cells[fi][0]  # cell C
        dfi = df[fi]
        efi = np.array(ef[fi])
        grad_uf[fi] = (uf[fi] - u[c1])/dfi*efi
        grad_vf[fi] = (vf[fi] - v[c1])/dfi*efi

    #######################
    ## pressure gradient ##
    #######################

    # face pressure
    for fi in range(n_faces):
        c1 = face_cells[fi][0]  # cell C
        c2 = face_cells[fi][1]  # cell F
        gfi = gf[fi]
        pf[fi] = gfi*p[c2] + (1 - gfi)*p[c1]
    for fi in outlet_faces:  # required
        pf[fi] = 0.0
    # cell pressure gradient
    for ci in range(n_cells):
        Vci = Vc[ci]
        grad_p[ci] = 0.0
        for j in range(3):
            fj = cell_faces[ci][j]
            Sfj = Sf[fj]
            nfj = nt[ci][j]*np.array(nf[fj])
            grad_p[ci] += (1/Vci)*pf[fj]*Sfj*nfj
    # face pressure gradient
    for fi in interior_faces:
        c1 = face_cells[fi][0]  # cell C
        c2 = face_cells[fi][1]  # cell F
        gfi = gf[fi]
        grad_pf_bar = gfi*grad_p[c2] + (1 - gfi)*grad_p[c1]

        dfi = df[fi]
        efi = np.array(ef[fi])
        grad_pf[fi] = grad_pf_bar + ((p[c2] - p[c1])/dfi - np.dot(grad_pf_bar, efi))*efi
    for fi in outlet_faces:
        c1 = face_cells[fi][0]  # cell C
        dfi = df[fi]
        efi = np.array(ef[fi])
        grad_pf[fi] = (pf[fi] - p[c1])/dfi*efi

    #######################
    ## momentum equation ##
    #######################

    for ci in interior_cells:
        aCu[ci] = 0.0
        bCu[ci] = 0.0
        aCv[ci] = 0.0
        bCv[ci] = 0.0
        for j in range(3):
            cj = cell_cells[ci][j]
            fj = cell_faces[ci][j]
            mfj = nt[ci][j]*mf[fj]
            Efj = Ef[fj]
            dfj = df[fj]
            efj = nt[ci][j]*np.array(ef[fj])
            Tfj = nt[ci][j]*np.array(Tf[fj])
            Sfj = Sf[fj]
            nfj = nt[ci][j]*np.array(nf[fj])

            aCu[ci] += max(mfj, 0) + nu*Efj/dfj
            aFu[ci, j] = -max(-mfj, 0) - nu*Efj/dfj
            if u[ci] == u[cj]:
                bCu[ci] += nu*np.dot(grad_uf[fj], Tfj) - pf[fj]*Sfj*nfj[0]
            else:
                uU = u[cj] - 2*np.dot(grad_u[ci], dfj*efj)
                uD = u[ci] + 2*np.dot(grad_u[cj], dfj*efj)
                rfp = (u[ci] - uU)/(u[cj] - u[ci])
                rfm = (u[cj] - uD)/(u[ci] - u[cj])
                psip = max(0, min(2*rfp, (rfp + 1)/2, 2))
                psim = max(0, min(2*rfm, (rfm + 1)/2, 2))
                bCu[ci] += nu*np.dot(grad_uf[fj], Tfj) - 0.5*(u[cj] - u[ci])*(psip*max(mfj, 0) + psim*max(-mfj, 0)) - pf[fj]*Sfj*nfj[0]

            aCv[ci] += max(mfj, 0) + nu*Efj/dfj
            aFv[ci, j] = -max(-mfj, 0) - nu*Efj/dfj
            if v[ci] == v[cj]:
                bCv[ci] += nu*np.dot(grad_vf[fj], Tfj) - pf[fj]*Sfj*nfj[1]
            else:
                vU = v[cj] - 2*np.dot(grad_v[ci], dfj*efj)
                vD = v[ci] + 2*np.dot(grad_v[cj], dfj*efj)
                rfp = (v[ci] - vU)/(v[cj] - v[ci])
                rfm = (v[cj] - vD)/(v[ci] - v[cj])
                psip = max(0, min(2*rfp, (rfp + 1)/2, 2))
                psim = max(0, min(2*rfm, (rfm + 1)/2, 2))
                bCv[ci] += nu*np.dot(grad_vf[fj], Tfj) - 0.5*(v[cj] - v[ci])*(psip*max(mfj, 0) + psim*max(-mfj, 0)) - pf[fj]*Sfj*nfj[1]
    for ci in boundary_cells:
        aCu[ci] = 0.0
        bCu[ci] = 0.0
        aCv[ci] = 0.0
        bCv[ci] = 0.0
        for j in range(3):
            cj = cell_cells[ci][j]
            fj = cell_faces[ci][j]
            mfj = nt[ci][j]*mf[fj]
            Efj = Ef[fj]
            dfj = df[fj]
            efj = nt[ci][j]*np.array(ef[fj])
            Tfj = nt[ci][j]*np.array(Tf[fj])
            Sfj = Sf[fj]
            nfj = nt[ci][j]*np.array(nf[fj])
            hfj = hf[fj]

            if fj in wall_faces + object_faces:
                aCu[ci] += (nu*Sfj/hfj)*(1 - nfj[0]**2)
                # aFu[ci, j] = 0.0
                bCu[ci] += (nu*Sfj/hfj)*(v[ci]*nfj[0]*nfj[1]) - pf[fj]*Sfj*nfj[0]

                aCv[ci] += (nu*Sfj/hfj)*(1 - nfj[1]**2)
                # aFv[ci, j] = 0.0
                bCv[ci] += (nu*Sfj/hfj)*(u[ci]*nfj[0]*nfj[1]) - pf[fj]*Sfj*nfj[1]
            elif fj in inlet_faces:
                aCu[ci] += nu*Efj/dfj
                # aFu[ci, j] = 0.0
                bCu[ci] += -(mfj - nu*Efj/dfj)*uf[fj] - pf[fj]*Sfj*nfj[0]

                aCv[ci] += nu*Efj/dfj
                # aFv[ci, j] = 0.0
                bCv[ci] += -(mfj - nu*Efj/dfj)*vf[fj] - pf[fj]*Sfj*nfj[1]
            elif fj in outlet_faces:
                aCu[ci] += mfj
                # aFu[ci, j] = 0.0
                # bCu[ci] += 0.0

                aCv[ci] += mfj
                # aFv[ci, j] = 0.0
                # bCv[ci] += 0.0
            else:
                aCu[ci] += max(mfj, 0) + nu*Efj/dfj
                aFu[ci, j] = -max(-mfj, 0) - nu*Efj/dfj
                if u[ci] == u[cj]:
                    bCu[ci] += nu*np.dot(grad_uf[fj], Tfj) - pf[fj]*Sfj*nfj[0]
                else:
                    uU = u[cj] - 2*np.dot(grad_u[ci], dfj*efj)
                    uD = u[ci] + 2*np.dot(grad_u[cj], dfj*efj)
                    rfp = (u[ci] - uU)/(u[cj] - u[ci])
                    rfm = (u[cj] - uD)/(u[ci] - u[cj])
                    psip = max(0, min(2*rfp, (rfp + 1)/2, 2))
                    psim = max(0, min(2*rfm, (rfm + 1)/2, 2))
                    bCu[ci] += nu*np.dot(grad_uf[fj], Tfj) - 0.5*(u[cj] - u[ci])*(psip*max(mfj, 0) + psim*max(-mfj, 0)) - pf[fj]*Sfj*nfj[0]

                aCv[ci] += max(mfj, 0) + nu*Efj/dfj
                aFv[ci, j] = -max(-mfj, 0) - nu*Efj/dfj
                if v[ci] == v[cj]:
                    bCv[ci] += nu*np.dot(grad_vf[fj], Tfj) - pf[fj]*Sfj*nfj[1]
                else:
                    vU = v[cj] - 2*np.dot(grad_v[ci], dfj*efj)
                    vD = v[ci] + 2*np.dot(grad_v[cj], dfj*efj)
                    rfp = (v[ci] - vU)/(v[cj] - v[ci])
                    rfm = (v[cj] - vD)/(v[ci] - v[cj])
                    psip = max(0, min(2*rfp, (rfp + 1)/2, 2))
                    psim = max(0, min(2*rfm, (rfm + 1)/2, 2))
                    bCv[ci] += nu*np.dot(grad_vf[fj], Tfj) - 0.5*(v[cj] - v[ci])*(psip*max(mfj, 0) + psim*max(-mfj, 0)) - pf[fj]*Sfj*nfj[1]

    for ci in range(n_cells):
        sum_u = 0.0
        sum_v = 0.0
        for j in range(3):
            cj = cell_cells[ci][j]
            sum_u += aFu[ci, j]*u[cj]
            sum_v += aFv[ci, j]*v[cj]
        us[ci] = (1/aCu[ci])*(-sum_u + bCu[ci])
        vs[ci] = (1/aCv[ci])*(-sum_v + bCv[ci])

    #######################
    ## pressure equation ##
    #######################

    for ci in range(n_cells):
        Vci = Vc[ci]
        DCu[ci] = Vci/aCu[ci]
        DCv[ci] = Vci/aCv[ci]
    for fi in range(n_faces):
        c1 = face_cells[fi][0]  # cell C
        c2 = face_cells[fi][1]  # cell F
        gfi = gf[fi]
        Dfu[fi] = gfi*DCu[c2] + (1 - gfi)*DCu[c1]
        Dfv[fi] = gfi*DCv[c2] + (1 - gfi)*DCv[c1]
        Dfi = np.array([Dfu[fi], Dfv[fi]])
        Sfi = Sf[fi]*np.array(nf[fi])
        Spfi = Dfi*Sfi
        efi = np.array(ef[fi])
        Epf[fi] = np.dot(Spfi, Spfi)/np.dot(Spfi, efi)

    for fi in interior_faces:
        c1 = face_cells[fi][0]  # cell C
        c2 = face_cells[fi][1]  # cell F
        gfi = gf[fi]
        usf_bar = gfi*us[c2] + (1 - gfi)*us[c1]
        vsf_bar = gfi*vs[c2] + (1 - gfi)*vs[c1]
        Vsf_bar = np.array([usf_bar, vsf_bar])
        grad_pf_bar = gfi*grad_p[c2] + (1 - gfi)*grad_p[c1]
        Dfi = np.array([Dfu[fi], Dfv[fi]])
        Sfi = Sf[fi]*np.array(nf[fi])
        Spfi = Dfi*Sfi
        msf[fi] = np.dot(Vsf_bar, Sfi) - np.dot(grad_pf[fi] - grad_pf_bar, Spfi)
    for fi in outlet_faces:
        c1 = face_cells[fi][0]  # cell C
        Vsfi = np.array([us[c1], vs[c1]])
        Sfi = Sf[fi]*np.array(nf[fi])
        msf[fi] = np.dot(Vsfi, Sfi)

    pp = np.zeros(n_cells)

    for ci in interior_cells:
        aCp[ci] = 0.0
        bCp[ci] = 0.0
        for j in range(3):
            fj = cell_faces[ci][j]
            msfj = nt[ci][j]*msf[fj]
            Epfj = Epf[fj]
            dfj = df[fj]

            aCp[ci] += -Epfj/dfj
            aFp[ci, j] = Epfj/dfj
            bCp[ci] += msfj
    for ci in boundary_cells:
        aCp[ci] = 0.0
        bCp[ci] = 0.0
        for j in range(3):
            fj = cell_faces[ci][j]
            msfj = nt[ci][j]*msf[fj]
            Epfj = Epf[fj]
            dfj = df[fj]

            # if fj in wall_faces + object_faces:
            #     aCp[ci] += 0.0
            #     aFp[ci, j] = 0.0
            #     bCp[ci] += 0.0
            # elif fj in inlet_faces:
            if fj in inlet_faces:
                # aCp[ci] += 0.0
                # aFp[ci, j] = 0.0
                bCp[ci] += msfj
            elif fj in outlet_faces:
                aCp[ci] += -Epfj/dfj
                # aFp[ci, j] = 0.0
                bCp[ci] += msfj
            else:
                aCp[ci] += -Epfj/dfj
                aFp[ci, j] = Epfj/dfj
                bCp[ci] += msfj

    for k in range(10):  # Gauss-Seidel
        for ci in range(n_cells):
            sum_p = 0.0
            for j in range(3):
                cj = cell_cells[ci][j]
                sum_p += aFp[ci, j]*pp[cj]
            pp[ci] = (1/aCp[ci])*(-sum_p + bCp[ci])

    ##################################
    ## pressure correction gradient ##
    ##################################

    # face pressure correction
    for fi in range(n_faces):
        c1 = face_cells[fi][0]  # cell C
        c2 = face_cells[fi][1]  # cell F
        gfi = gf[fi]
        ppf[fi] = gfi*pp[c2] + (1 - gfi)*pp[c1]
    for fi in outlet_faces:  # required
        ppf[fi] = 0.0
    # cell pressure correction gradient
    for ci in range(n_cells):
        Vci = Vc[ci]
        grad_pp[ci] = 0.0
        for j in range(3):
            fj = cell_faces[ci][j]
            Sfj = Sf[fj]
            nfj = nt[ci][j]*np.array(nf[fj])
            grad_pp[ci] += (1/Vci)*ppf[fj]*Sfj*nfj

    ######################
    ## update variables ##
    ######################

    # velocity and pressure
    for ci in range(n_cells):
        DCi = np.array([DCu[ci], DCv[ci]])
        Vsi = np.array([us[ci], vs[ci]])
        Vi = Vsi - DCi*grad_pp[ci]
        u[ci] = Vi[0]
        v[ci] = Vi[1]
        p[ci] = p[ci] + prf*pp[ci]
    # mass flow rate
    for fi in interior_faces:
        c1 = face_cells[fi][0]  # cell C
        c2 = face_cells[fi][1]  # cell F
        Epfi = Epf[fi]
        dfi = df[fi]
        mf[fi] = msf[fi] - Epfi*(pp[c2] - pp[c1])/dfi
    for fi in outlet_faces:
        c1 = face_cells[fi][0]  # cell C
        Epfi = Epf[fi]
        dfi = df[fi]
        mf[fi] = msf[fi] - Epfi*(ppf[fi] - pp[c1])/dfi

    ######################
    ## check continuity ##
    ######################

    mf_out = 0.0
    for fi in outlet_faces:
        mf_out += mf[fi]
    err[n] = (mf_out - mf_in)/mf_in * 100
    print('Iteration = %d, Global Mass Imbalance (%%) = %.5f' % (n, err[n]))

###############
## plot data ##
###############

# plot global mass imbalance
plt.plot(err)
plt.title('Global Mass Imbalance', fontweight='bold')
plt.xlabel('Iteration')
plt.ylabel('Imbalance (%)')
plt.grid('on')
plt.show()

# plot speed
np.save('u', u)
np.save('v', v)
np.save('p', p)
V = np.sqrt(u*u + v*v)

plt.figure(figsize=(10, 5))
cmap = plt.get_cmap('jet')
norm = clr.Normalize(vmin=min(V), vmax=max(V))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
plt.colorbar(sm, orientation='horizontal', label='Speed [m/s]', shrink=0.5)

for i in range(n_cells):
    n1 = nodes[cells[i][0]]  # cell node 1
    n2 = nodes[cells[i][1]]  # cell node 2
    n3 = nodes[cells[i][2]]  # cell node 3
    plt.fill([n1[0], n2[0], n3[0], n1[0]], [n1[1], n2[1], n3[1], n1[1]], c=cmap(norm(V[i])))

for i in boundary_faces:
    n1x = nodes[faces[i][0]][0]  # face node 1 x-coordinate
    n1y = nodes[faces[i][0]][1]  # face node 1 y-coordinate
    n2x = nodes[faces[i][1]][0]  # face node 2 x-coordinate
    n2y = nodes[faces[i][1]][1]  # face node 2 y-coordinate
    plt.plot([n1x, n2x], [n1y, n2y], c='black', lw=1.0)

plt.title('Re = %d, U = %.1f m/s, ν = %.3f m²/s' % (U*L/nu, U, nu), fontweight='bold')
plt.xlabel('Channel Length [m]')
plt.ylabel('Channel Height [m]')
plt.axis('equal')
plt.show()
