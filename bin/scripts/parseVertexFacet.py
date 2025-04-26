#!/usr/bin/env python3

###########
# parseVertexFacet.py -- SRS 2/19/16
# =========
# Parses vertex-facet shape model file and outputs pkdgrav wall file.
# Specific to input file format (x,y,z in km):
#     ["v" x y z] .... ["f" v1 v2 v3] ....
# Note: to use output, MAX_NUM_WALLS in walls.h may need to be increased.
###########

from os.path import isfile
import sys


def syntax():
    sys.stderr.write('syntax: ' + sys.argv[0] +
                     ' [shape file (vertex-facet)] [walls data output file]\n')
    sys.exit(1)


# Syntax checks

err = 0
if len(sys.argv) != 3:
    syntax()
if not isfile(sys.argv[1]):
    sys.stderr.write('error: unable to find specified shape file "' +
                     sys.argv[1] + '\n')
    err = 1
if isfile(sys.argv[2]):
    sys.stderr.write('error: specified walls data output file "' +
                     sys.argv[2] + '" already exists\n')
    err = 1
if (err):
    syntax()

# Parse file, write output

wallsDAT = open(sys.argv[2], 'w')
wallsDAT.write('lengthunit 6.6845871226705983e-9 # au in km\n')
wallsDAT.write('defaults epsn 0.0 color 170\n')

xV = ["ERROR"]
yV = ["ERROR"]
zV = ["ERROR"]
Fa = []
Fb = []
Fc = []
with open(sys.argv[1], 'r') as vfDAT:
    for line in vfDAT:
        row = (line.split(" "))
        c = row[0]
        if c == "v":
            xV.append(float(row[1]))
            yV.append(float(row[2]))
            zV.append(float(row[3]))
        elif c == "f":
            Fa = int(row[1])
            Fb = int(row[2])
            Fc = int(row[3])
            vV1 = [xV[Fa], yV[Fa], zV[Fa]]
            vV2 = [xV[Fb], yV[Fb], zV[Fb]]
            vV3 = [xV[Fc], yV[Fc], zV[Fc]]
            wallsDAT.write('wall type triangle\n')
            wallsDAT.write('origin %g %g %g\n' % (vV1[0], vV1[1], vV1[2]))
            wallsDAT.write('vertex1 %g %g %g\n' %
                           (vV2[0] - vV1[0], vV2[1] - vV1[1], vV2[2] - vV1[2]))
            wallsDAT.write('vertex2 %g %g %e\n' %
                           (vV3[0] - vV1[0], vV3[1] - vV1[1], vV3[2] - vV1[2]))
        else:
            break
wallsDAT.close()
