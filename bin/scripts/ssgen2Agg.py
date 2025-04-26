#!/usr/bin/env python3
import numpy as np
import sys
import subprocess
import argparse


def vectorDot(v1, v2):
    return (v1[0]*v2[0])+(v1[1]*v2[1])+(v1[2]*v2[2])


# Rotate a point (vPoint) about an axis (vAxis) by a an angle (dAngle)
def vectorRotate(vPoint, vAxis, dAngle):
    sa = np.sin(dAngle)
    ca = np.cos(dAngle)
    dot = vectorDot(vPoint, vAxis)
    x = vPoint[0]
    y = vPoint[1]
    z = vPoint[2]
    u = vAxis[0]
    v = vAxis[1]
    w = vAxis[2]

    vPoint[0] = u*dot+(x*(v*v+w*w) - u*(v*y+w*z))*ca + (-w*y+v*z)*sa
    vPoint[1] = v*dot+(y*(u*u+w*w) - v*(u*x+w*z))*ca + (w*x-u*z)*sa
    vPoint[2] = w*dot+(z*(u*u+v*v) - w*(u*x+v*y))*ca + (-v*x+u*y)*sa

    return vPoint


# Rotate a rigid body (assumed to be centered at the origin)
# along the x and y by random angles
def rotateRigid(body):
    xAxis = [1, 0, 0]
    yAxis = [0, 1, 0]
    dAngleX = 2*np.pi*np.random.rand()
    dAngleY = 2*np.pi*np.random.rand()
    for i in np.arange(body.shape[0]):
        body[i, 4:7] = vectorRotate(body[i, 4:7], xAxis, dAngleX)
        body[i, 4:7] = vectorRotate(body[i, 4:7], yAxis, dAngleY)
    return body


# Generate a Dumbbell aggregate within a Bounding Volume
def genDB(bodyBoundVol, sep, radius, index, neg_index, color):
    lengthunit = 6.68458712e-14
    refPoint = bodyBoundVol[0, 4:7]
    DB = bodyBoundVol
    num_particles = 2
    for i in np.arange(num_particles-1):
        DB = np.concatenate((DB, bodyBoundVol))
    r = radius*lengthunit
    s = sep*lengthunit
    DB[0, 4:7] = [0, 0, -s/2.0]
    DB[1, 4:7] = [0, 0, +s/2.0]
    DB[:, 0] = np.arange(DB.shape[0])+index
    DB[:, 1] = neg_index
    DB[:, 2] = bodyBoundVol[0, 2]*r*r*r / \
        (bodyBoundVol[0, 3]*bodyBoundVol[0, 3]*bodyBoundVol[0, 3])
    DB[:, 3] = r
    DB = rotateRigid(DB)
    DB[:, 4:7] = DB[:, 4:7]+refPoint
    if color == True:
        DB[:, 13] = 2
    elif color == False:
        DB[:, 13] = 1
    return DB


# Generate a Tetrahedron aggregate within a Bounding Volume
def genTetrahedron(bodyBoundVol, sep, radius, index, neg_index, color):
    lengthunit = 6.68458712e-14
    refPoint = bodyBoundVol[0, 4:7]
    tetra = bodyBoundVol
    num_particles = 4
    for i in np.arange(num_particles-1):
        tetra = np.concatenate((tetra, bodyBoundVol))
    r = radius*lengthunit
    s = sep*lengthunit
    tetra[0, 4:7] = [s/2.0, 0, -s/(2.0*np.sqrt(2))]
    tetra[1, 4:7] = [-s/2.0, 0, -s/(2.0*np.sqrt(2))]
    tetra[2, 4:7] = [0, -s/2.0, s/(2.0*np.sqrt(2))]
    tetra[3, 4:7] = [0, s/2.0, s/(2.0*np.sqrt(2))]
    tetra[:, 0] = np.arange(tetra.shape[0])+index
    tetra[:, 1] = neg_index
    tetra[:, 2] = bodyBoundVol[0, 2]*r*r*r / \
        (bodyBoundVol[0, 3]*bodyBoundVol[0, 3]*bodyBoundVol[0, 3])
    tetra[:, 3] = r
    tetra = rotateRigid(tetra)
    tetra[:, 4:7] = tetra[:, 4:7]+refPoint
    if color == True:
        tetra[:, 13] = 3
    elif color == False:
        tetra[:, 13] = 1
    return tetra


# Generate a 4 particle Rod aggregate within a Bounding Volume
def gen4Rod(bodyBoundVol, sep, radius, index, neg_index, color):
    lengthunit = 6.68458712e-14
    refPoint = bodyBoundVol[0, 4:7]
    rod = bodyBoundVol
    num_particles = 4
    for i in np.arange(num_particles-1):
        rod = np.concatenate((rod, bodyBoundVol))
    r = radius*lengthunit
    s = sep*lengthunit
    rod[0, 4:7] = [0, 0, -3*s/2.0]
    rod[1, 4:7] = [0, 0, -s/2.0]
    rod[2, 4:7] = [0, 0, +s/2.0]
    rod[3, 4:7] = [0, 0, +3*s/2.0]
    rod[:, 0] = np.arange(rod.shape[0])+index
    rod[:, 1] = neg_index
    rod[:, 2] = bodyBoundVol[0, 2]*r*r*r / \
        (bodyBoundVol[0, 3]*bodyBoundVol[0, 3]*bodyBoundVol[0, 3])
    rod[:, 3] = r
    rod = rotateRigid(rod)
    rod[:, 4:7] = rod[:, 4:7]+refPoint
    if color == True:
        rod[:, 13] = 5
    elif color == False:
        rod[:, 13] = 1
    return rod


# Generate a Diamond shape aggregate within a Bounding Volume
def genDiamond(bodyBoundVol, sep, radius, index, neg_index, color):
    lengthunit = 6.68458712e-14
    refPoint = bodyBoundVol[0, 4:7]
    diamond = bodyBoundVol
    num_particles = 4
    for i in np.arange(num_particles-1):
        diamond = np.concatenate((diamond, bodyBoundVol))
    r = radius*lengthunit
    s = sep*lengthunit
    diamond[0, 4:7] = [0, 0, -s]
    diamond[1, 4:7] = [0, 0, +s]
    diamond[2, 4:7] = [-s/2., 0, 0]
    diamond[3, 4:7] = [+s/2., 0, 0]
    diamond[:, 0] = np.arange(diamond.shape[0])+index
    diamond[:, 1] = neg_index
    diamond[:, 2] = bodyBoundVol[0, 2]*r*r*r / \
        (bodyBoundVol[0, 3]*bodyBoundVol[0, 3]*bodyBoundVol[0, 3])
    diamond[:, 3] = r
    diamond = rotateRigid(diamond)
    diamond[:, 4:7] = diamond[:, 4:7]+refPoint
    if color == True:
        diamond[:, 13] = 12
    elif color == False:
        diamond[:, 13] = 1
    return diamond


# Generate a Cube shape aggregate within a Bounding Volume
def genCube(bodyBoundVol, sep, radius, index, neg_index, color):
    lengthunit = 6.68458712e-14
    refPoint = bodyBoundVol[0, 4:7]
    cube = bodyBoundVol
    num_particles = 8
    for i in np.arange(num_particles-1):
        cube = np.concatenate((cube, bodyBoundVol))
    r = radius*lengthunit
    s = sep*lengthunit
    cube[0, 4:7] = [-s/2.0, -s/2.0, -s/2.0]
    cube[1, 4:7] = [s/2.0, -s/2.0, -s/2.0]
    cube[2, 4:7] = [-s/2.0, s/2.0, -s/2.0]
    cube[3, 4:7] = [s/2.0, s/2.0, -s/2.0]
    cube[4, 4:7] = [-s/2.0, -s/2.0, s/2.0]
    cube[5, 4:7] = [s/2.0, -s/2.0, s/2.0]
    cube[6, 4:7] = [-s/2.0, s/2.0, s/2.0]
    cube[7, 4:7] = [s/2.0, s/2.0, s/2.0]
    cube[:, 0] = np.arange(cube.shape[0])+index
    cube[:, 1] = neg_index
    cube[:, 2] = bodyBoundVol[0, 2]*r*r*r / \
        (bodyBoundVol[0, 3]*bodyBoundVol[0, 3]*bodyBoundVol[0, 3])
    cube[:, 3] = r
    cube = rotateRigid(cube)
    cube[:, 4:7] = cube[:, 4:7]+refPoint
    if color == True:
        cube[:, 13] = 7
    elif color == False:
        cube[:, 13] = 1
    return cube

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                        description = 'Generates a field of \
                        aggregates of a given shape and writes to ssgen.ss.',
                        usage = 'ssgen2Agg.py shape-code radius density '
                        'lx ly lz N_aggs [-s sep] [-i idx] [-c]  ')
    parser.add_argument('shape', metavar='shape-code', type=int,
                        help='0 = Tetrahedron, 1 = 4 particle Rod, 2 = Dumbbell,\
                        3 = Diamond, 4 = Cube, 5 = Mix')
    parser.add_argument('radius', metavar='radius', type=float,
                        help='Size of individual particles in the aggregate (cm)')
    parser.add_argument('density', metavar='density', type=float,
                        help='Individual particle density (g/cc)')
    parser.add_argument('region', metavar=('lx ly lz'), type=float, nargs=3,
                        help='Simulation region dimensions (m)')
    parser.add_argument('N_aggs', type=int, help='The number \
                         of aggregates, not particles.')
    parser.add_argument('-s', metavar='sep', type=float,
                        help='Separation between particles in the aggregate,\
                        in multiples of the particle radius (default sqrt(3)).')
    parser.add_argument('-i', metavar='idx', type=int,
                        help='Give starting index for aggs (default 0)')
    parser.add_argument('-c', action='store_true', help='Color aggregates by '
                        'their shapes (off by default).')
    # Some optional arguments from ssgen. JCM 10/1/21
    group = parser.add_mutually_exclusive_group()
    parser.add_argument('-n', metavar='max', type=int, help='The maximum number \
                        of passes allowed (default 500).')
    group.add_argument('-e', action='store_true', help='Ellipsoidal, not \
                        rectangular, region.')
    group.add_argument('-y', action='store_true', help='Cylindrical region \
                        (lx,ly are cross-section dimensions).')

    args = parser.parse_args()

    # Run ssgen with user-input parameters
    shape = args.shape
    radius = args.radius
    density = args.density
    lx = args.region[0]
    ly = args.region[1]
    lz = args.region[2]
    N_aggs = args.N_aggs
    if args.i is not None:
        idx = args.idx
    else: 
        idx = 0
    if args.s is not None:
        sep_coeff = args.s
    else:
        sep_coeff = np.sqrt(3)
    color = args.c

    # Adjust aggregate "radius" for ssgen depending on shape of particle
    mix = 0
    if shape == 0:
        agg_radius = (1.0 + 2.0 / np.sqrt(3.0))*radius
    elif shape == 1:
        agg_radius = 4*radius
    elif shape == 2:
        agg_radius = 2*radius
    elif shape == 3:
        agg_radius = 2.5*radius
    elif shape == 4:
        agg_radius = (1.0 + np.sqrt(3.0))*radius
    elif shape == 5:
        agg_radius = 4*radius
        mix = 1
    else:
        print("Unavailable shape code")
        print("shape codes:")
        print("0 = Tetrahedron")
        print("1 = 4 particle Rod")
        print("2 = Dumbbell")
        print("3 = Diamond")
        print("4 = Cube")
        print("5 = Mix")
        sys.exit(1)

    # Run ssgen and convert to bt file
    ssgenCommand = "ssgen -r "+str(agg_radius)+" -d "+str(density)+" -x " + \
        str(lx)+" -y "+str(ly)+" -z "+str(lz)+" "+str(N_aggs)   
    # Optional arguments to pass to ssgen
    if args.e is True:
        ssgenCommand += " -e"
    elif args.y is True:
        ssgenCommand += " -c"
    if args.n is not None:
        ssgenCommand += " -n "+str(args.n)
    
    ss2btCommand = "ss2bt ssgen.ss"
    process = subprocess.call(ssgenCommand.split())
    if process == 1:
        print("ERROR: SSGEN failed.")
        sys.exit(1)
    process = subprocess.call(ss2btCommand.split())
    if process == 1:
        print("ERROR: Could not convert ssgen.ss to ssgen.bt")
        sys.exit(1)

    # negative index, count with base 3
    L = -1-idx
    #
    a = 0

    particles = np.genfromtxt("ssgen.bt")
    particles = np.reshape(particles, (N_aggs, 14))
    f1 = open('ssgen.bt', 'w')

    sep = sep_coeff*radius  # in cm
    for i in range(N_aggs):
        bodyBoundVol = particles[i, :]
        bodyBoundVol = np.reshape(bodyBoundVol, (1, 14))
        if mix == 1:
            shape = np.random.randint(0,5)

        if shape == 0:
            body = genTetrahedron(bodyBoundVol, sep, radius, a, L, color)
            L = L-1
            a = a+4
        if shape == 1:
            body = gen4Rod(bodyBoundVol, sep, radius, a, L, color)
            L = L-1
            a = a+4
        if shape == 2:
            body = genDB(bodyBoundVol, sep, radius, a, L, color)
            L = L-1
            a = a+2
        if shape == 3:
            body = genDiamond(bodyBoundVol, sep, radius, a, L, color)
            L = L-1
            a = a+4
        if shape == 4:
            body = genCube(bodyBoundVol, sep, radius, a, L, color)
            L = L-1
            a = a+8
        for j in np.arange(body.shape[0]):
            f1.write("%d " % body[j, 0]+"%d " % body[j, 1] +
                     "%.16e " % body[j, 2] +
                     "%.16e " % body[j, 3]+"%.16e " % body[j, 4] +
                     "%.16e " % body[j, 5]+"%.16e " % body[j, 6] +
                     "%.16e " % body[j, 7]+"%.16e " % body[j, 8] +
                     "%.16e " % body[j, 9]+"%.16e " % body[j, 10] +
                     "%.16e " % body[j, 11]+"%.16e " % body[j, 12] +
                     "%d" % body[j, 13]+"\n")

    f1.close()
    bt2ssCommand = "bt2ss ssgen.bt"
    process = subprocess.call(bt2ssCommand.split())
    if process == 1:
        print("ERROR: Could not convert ssgen.bt to ssgen.ss")
        sys.exit(1)
    else:
        print("Success: Created Aggregate")
