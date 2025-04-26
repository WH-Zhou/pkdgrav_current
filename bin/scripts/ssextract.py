#!/usr/bin/env python3

import argparse
import numpy as np
import os
import sys

try:
    import ssio
except ModuleNotFoundError:
    print("Error: ssio library not found.")
    print("Add location of ssio.py to PYTHONPATH environment variable.")
    print("Or copy ssio.py to this directory.")
    print("Look for ssio.py in the source code bin/scripts directory.")
    sys.exit(1)


# Define functions to start/stop print commands from displaying
def stop_print():
    sys.stdout = open(os.devnull, 'w')


def start_print():
    sys.stdout = sys.__stdout__


# Define function to recenter the object wrt the center of mass
def cm(x, y, z, m):
    M = np.sum(m)
    x_cm = np.sum(m*x/M)
    y_cm = np.sum(m*y/M)
    z_cm = np.sum(m*z/M)
    return x_cm, y_cm, z_cm


# Parse Arguments
parser = argparse.ArgumentParser(
        description='Extracts particles of a certain color or in a specified '
        'region and rewrites the .ss files in the specified directory. '
        'NOTE: Particle radius is not taken into account when extracting '
        'particles from a region.',
        usage='ssextract.py [-h] [-q] [-o] [-u units] [-f] ([-c color] '
        '[-r x y z | -e x y z]) path ss_files',
        epilog='NOTE: The user must use at least one of -c, -e, or -r; '
        '-c can be used in conjunction with -e and -r, but -e and -r are '
        'mutually exclusive.\n\nColors: 0=BLACK 1=WHITE 2=RED 3=GREEN '
        '4=BLUE 5=YELLOW 6=MAGENTA 7=CYAN\n8=GOLD 9=PINK 10=ORANGE '
        '11=KHAKI 12=VIOLET 13=MAROON 14=AQUA 15=NAVY\n16=BLACK '
        '17-254=INCREASING GRAY 255=WHITE -N=STUCK TO WALL 1-N '
        '>255=24-BIT\n')
parser.add_argument("-c", "--color", type=int, nargs=1, help='An integer'
                    ' denoting the desired color of particles to be extracted'
                    ' (color-integer relationship detailed below). Must be '
                    'used if -r or -e not used.')
shapes = parser.add_mutually_exclusive_group()
shapes.add_argument("-r", "--rectangle", metavar=('x', 'y', 'z'),  type=float,
                    nargs=3, help="Full lengths of the sides of the "
                    "rectangular shape to be extracted, centered on the "
                    "origin. Must be used if -e or -c not used.")
shapes.add_argument("-e", "--ellipsoid", metavar=('x', 'y', 'z'), type=float,
                    nargs=3, help="Semi-axis lengths of the ellipsoidal shape "
                    "to be extracted (radii, not diameters), measured from "
                    "the origin. Must be used if -r or -c not used.")
parser.add_argument("path", help="Path of desired output specified relative "
                    "to root or current working directory. Must be empty "
                    "unless -f option used.")
parser.add_argument("ss_files", nargs='+', help='List of .ss files to be '
                    'worked on.')
parser.add_argument("-o", "--origin", action="store_false", default=True,
                    help="To be used with -e or -r options. Give bounds "
                    "relative to center of mass and recenter object to the "
                    "origin. If this option is not used, bounds must be "
                    "specified with respect to the origin.")
parser.add_argument("-a", "--rem_aggs", action="store_true", default=False,
                    help="Remove all aggregate constituents if any member "
                    "crosses specified boundaries. By default, only the "
                    "particle crossing a boundary specified with -e or -r will "
                    "be removed from the agg.")
parser.add_argument("-u", "--units", nargs='?', default='km', help="To be used"
                    " with -e or -r options: change units for input lengths. "
                    "Default = Kilometers, Meters = 'm', AU = 'AU'.")
parser.add_argument("-q", "--quiet", action="store_true", help="Silence "
                    "command line outputs.")
parser.add_argument("-f", "--force", action="store_true", help="*USE WITH "
                    "CAUTION* Remove files in PATH and force outputs to be "
                    "written into existing directory specified by PATH.")
args = parser.parse_args()

if args.color is None and args.rectangle is None and args.ellipsoid is None:
    print('usage: extract.py [-h] [-q] [-f] ([-c color] [-r x y z | -e x y z])'
          ' path ss_files\nssextract.py: error: Must input geometry and/or '
          'color.')
    sys.exit(1)

if args.units != 'km' and args.units != 'AU' and args.units != 'm':
    print('usage: extract.py [-h] [-q] [-f] ([-c color] [-r x y z | -e x y z])'
          ' path ss_files\nssextract.py: error: Input units as either km, m, '
          'or AU (Default = km).')
    sys.exit(1)

# Define variables based on command line arguments
path = args.path
km = 1.496e8		# km per AU
me = 1.496e11		# m per AU

if args.ellipsoid is not None:
    a = args.ellipsoid[0]
    b = args.ellipsoid[1]
    c = args.ellipsoid[2]
    if args.units == 'km':
        a = a/km
        b = b/km
        c = c/km
    elif args.units == 'm':
        a = a/me
        b = b/me
        c = c/me

if args.rectangle is not None:
    dx = args.rectangle[0]/2
    dy = args.rectangle[1]/2
    dz = args.rectangle[2]/2
    if args.units == 'km':
        dx = dx/km
        dy = dy/km
        dz = dz/km
    elif args.units == 'm':
        dx = dx/me
        dy = dy/me
        dz = dz/me

# Check to see if directory exists and is empty
if os.path.isdir(path):
    if args.force:					# Force writing to non-empty directory
        print('Removing existing files in %s' % path)
        for file in os.listdir(path):
            if os.path.isdir(file):
                print('Error: Remove all subdirectories from %s/' % path)
                sys.exit(2)
            else:
                os.remove(os.path.join(path, file))
    else:
        if len(os.listdir(path)) != 0:  # Exit if directory not empty
            print('Error: Target directory not empty.')
            sys.exit(3)
else:  # Create desired directory if nonexistent
    print ('Directory does not exist. Creating directory %s' % path)
    os.mkdir(path)

# Extract specified particles
for file in args.ss_files:
    try:
        header, data = ssio.read_SS(file, 'y')		# Read ss-file
    except ssio.ReadError as msg:		# Check for read errors (ssio_py.pdf)
        print('Read error:', msg)
        sys.exit(4)
    m = data[2, :]				# Read in mass and position data
    x = data[4, :]
    y = data[5, :]
    z = data[6, :]
    rem = np.empty(0, dtype = 'int64')		# Initialize array of columns to delete
    if args.color:
        for col in range(0, header[1]):
            # Find particles with color different from input
            if data[13, col] != args.color[0]:
                # Create an array of particles (columns) to delete
                rem = np.append(rem, col)
    if args.ellipsoid or args.rectangle:
        if args.origin is False:			# Change to Center of Mass coords
            x_cm, y_cm, z_cm = cm(x, y, z, m)
            x = x-x_cm
            y = y-y_cm
            z = z-z_cm
    if args.ellipsoid:
        for col in range(0, header[1]):
            # Find particles outside ellipsoidal region given by a, b, c
            if (x[col]/a)**2 + (y[col]/b)**2 + (z[col]/c)**2 > 1:
                rem = np.append(rem, col)
    elif args.rectangle:
        for col in range(0, header[1]):
            # Find particles outside of rectangular region given by x, y, z
            if np.absolute(x[col]) > dx or np.absolute(y[col]) > dy or \
                    np.absolute(z[col]) > dz:
                rem = np.append(rem, col)
    if args.rem_aggs:
        for particle in range(len(rem)):
            # Pick particles w/ original index matching any already found
            rem = np.append(rem, data[0, data[1] == data[1,rem[particle]]])
        rem = np.unique(rem)	# Sorts & solves multiple instances of same particle in rem
    data = np.delete(data, rem, 1)

# Rewrite ss-files to new directory
    if args.quiet:						# Silence outputs
        stop_print()
    try:
        # Rewrite ss-file to new directory
        ssio.write_SS(data, path+'/'+file, time=header[0])
    except ssio.WriteError as msg:
        print('Write error:', msg)
        sys.exit(5)
    if args.quiet:
        start_print()					# Unsilence outputs
    elif not args.quiet:
        print('Extracting file: '+file)
