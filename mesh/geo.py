import numpy as np

L = 10.         # domain length
H = 3.          # domain height
xle = 1.5       # leading edge x-coordinate
yle = 1.5       # leading edge y-coordinate
esf = .1        # mesh size factor

XY = np.loadtxt('mesh/body.txt')
X, Y = XY[:, 0] + xle, XY[:, 1] + yle

geo = open('mesh/mesh.geo', 'w')

geo.write('\n// domain points\n\n') # counter-clockwise
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (1, 0., 0.))
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (2, L, 0.))
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (3, L, H))
geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (4, 0., H))

geo.write('\n// body points\n\n') # clockwise
for i in range(len(X)):
    geo.write('Point(%d) = {%.10f, %.10f, 0, 1.0};\n' % (i + 5, X[i], Y[i]))

geo.write('\n// domain lines\n\n')
geo.write('Line(1) = {1, 2};\n')  # bottom
geo.write('Line(2) = {2, 3};\n')  # right
geo.write('Line(3) = {3, 4};\n')  # top
geo.write('Line(4) = {4, 1};\n')  # left

geo.write('\n// body lines\n\n')
for i in range(len(X) - 1):
    geo.write('Line(%d) = {%d, %d};\n' % (i + 5, i + 5, i + 6))
geo.write('Line(%d) = {%d, %d};\n' % (len(X) + 4, len(X) + 4, 5))

geo.write('\n// curve loops\n\n')
geo.write('Curve Loop(1) = {1, 2, 3, 4};\n')
body = ', '.join(map(str, list(range(5, len(X) + 5))))
geo.write('Curve Loop(2) = {' + body + '};\n')

geo.write('\n// plane surface\n\n') # domain boundary before body boundary
geo.write('Plane Surface(1) = {1, 2};\n')

geo.write('\n// physical groups\n\n')
geo.write('Physical Curve("OUTLET", %d) = {2};\n' % (len(X) + 5))
geo.write('Physical Curve("INLET", %d) = {4};\n' % (len(X) + 6))
geo.write('Physical Curve("WALL", %d) = {1, 3, ' % (len(X) + 7) + body + '};\n')

geo.write('\n// size field\n\n')
geo.write('Mesh.MeshSizeFactor = %.10f;\n' % esf)
