#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Mon Jun  8 10:16:20 2020

@author: dcr

Utility to make a version of the pkdgrav code base that omits certain macros.
Designed for distributing code beyond the project team.
'''

import argparse
import glob
import os
import shutil
import subprocess
import sys


def myignore(dir, contents):
    # custom utility function for shutil.copytree() to skip certain items...
    ignore = []
    for i in [args.export, '.git', '.gitignore', '__pycache__', 'build.log',
              'export_strip.py', 'export_strip_py.tex', 'export_strip_py.pdf']:
        if i in contents:
            ignore.append(i)
    return ignore


# macros that are removed automatically...
# (PHOBOS, REPEL_MARK_II, SIMPLE_GAS_DRAG, SPECIAL_PARTICLES ignored for now)

auto = ['NUMREC']  # these generate informative messages
extras = [  # these do not
    'CHARGE', 'DEM_TIDAL_LOCAL', 'DEM_TIDAL_SPACE', 'GLASS_BEADS',
    '_GLASS_BEADS', 'GR_DRAG', 'JOHNNY', 'RORY_EXTENSION', 'ROT_FRAME'
    ]

# macro dependencies to assist user... (using dict for convenient access)

depend = dict()
depend['COLLISIONS'] = {
    'WALLS', 'AGGS', 'SPRINGS', 'SLIDING_PATCH', 'DEM', 'SPINUP',
    'SEISMIC_WAVE', 'RUBBLE_ZML', 'COLLMOD_ZML', 'SPECIAL_PARTICLES'
    }
depend['WALLS'] = {'DEM_WALLS_REACT'}
depend['AGGS'] = {'FAST_AGGS'}
depend['DEM'] = {
    'DEM_FIXED_BALL', 'DEM_ROTATION_DASHPOT', 'DEM_TRACK_ORIENT',
    'DEM_COHESION', 'DEM_TWOLAYERS', 'DEM_WALLS_REACT', 'DEM_OVERLAP_OUTPUT',
    'DEM_DIAG'
    }
depend['DEM_COHESION'] = {'DEM_ROTATION_DASHPOT'}
depend['SPINUP'] = {'AGGS'}
depend['SEISMIC_WAVE'] = {'WALLS'}

# get user argument(s)...

parser = argparse.ArgumentParser(description='remove pkdgrav macro blocks',
                                 epilog='NUMREC is removed automatically.')
parser.add_argument('macro', nargs='*', help='macro(s) to remove')
parser.add_argument('-e', '--export', default='export',
                    help='name of export directory (default "%(default)s")')
parser.add_argument('-v', '--verbose', action='store_true', help='be verbose')
args = parser.parse_args()

# construct final list of user-specified macros to remove...

macros = []
# args.macro can be NoneType (0 items), str (1 item), or list (> 1 item)...
if type(args.macro) == str:
    macros.append(args.macro)
elif type(args.macro) == list:
    macros.extend(args.macro)
for m in auto:
    if m not in macros:
        macros.append(m)  # add automatic macros that have messaging

# check macro list...

if not os.path.isfile('Makefile.in'):
    print('Makefile.in not found in current directory.\n'
          'Be sure to run this script at the top of the codebase.')
    sys.exit(1)

try:
    with open('Makefile.in') as f:
        # read file as list split on whitespace, then make set by splitting
        # on "=" (keeping the leftmost item) and removing any lead "#"...
        # (this will isolate items like "USE_COLLISIONS", etc.)
        t = set([x.split('=')[0].lstrip('#') for x in f.read().split()])
except IOError:
    print('Unable to access Makefile.in.')
    sys.exit(1)
bCapsMsgGiven = False
for m in list(macros):  # use list() here to make new copy so remove() works
    if 'USE_' + m not in t:
        print(f'USE_{m} not found in Makefile.in...skipping...')
        macros.remove(m)
        if not bCapsMsgGiven and not m.isupper():
            print('   (note most macros are given in ALL CAPS)')
            bCapsMsgGiven = True

# apply dependencies...

for m in macros:
    if m in depend:
        macros.extend(depend[m])

macros = sorted(list(set(macros)))  # to remove duplicates (and reorder)

# finally, some special derived dependencies (no messaging)...

if 'AGGS' in macros or 'SLIDING_PATCH' in macros:
    extras.append('AGGS_IN_PATCH')
if 'COLLISIONS' in macros:
    extras.append('SSIO_USE_MPI')
if 'RUBBLE_ZML' in macros or 'COLLMOD_ZML' in macros:
    extras.append('ORIGIN_HISTOGRAM')

# prepare export...

top = os.path.basename(os.getcwd())
print(f'This script will run "make clean" in "{top}".')
input('Press return to proceed...')

result = subprocess.run('make clean', capture_output=True, shell=True)
if result.returncode != 0:
    print(f'Problem cleaning "{top}":\n{result.stdout.decode("utf-8")}')
    sys.exit(1)

try:
    os.mkdir(args.export)
except FileExistsError:
    print(f'"{args.export}" already exists -- remove or rename to continue.')
    sys.exit(1)

if args.verbose:
    print('Copying code base...', flush=True)

# note that "symlinks=True" is not needed because unifdef will overwrite
# symlinks anyway (so a few files will be duplicated in the export)...
shutil.copytree('.', f'{args.export}/{top}', ignore=myignore)

if args.verbose:
    print('Macros to remove:', end='', flush=True)
    for m in macros + extras:
        print(f' {m}', end='', flush=True)
    print()

# use "unifdef" tool to remove unwanted macros from .c & .h files...

files = glob.glob(f'{args.export}/{top}/src/pkdgrav/*.c')
files.extend(glob.glob(f'{args.export}/{top}/src/pkdgrav/*.h'))

strip = ''
for m in macros + extras:
    strip += f' -U{m}'
for f in files:
    if args.verbose:
        print(f'Processing {f}...', flush=True)
    result = subprocess.run(f'unifdef{strip} -B -o "{f}" "{f}"',
                            capture_output=True, shell=True)
    if result.returncode == 2:
        print(f'Error processing {f}:\n{result.stdout.decode("utf-8")}')
        sys.exit(1)

# final informative messages...

print(f'\nTo finish, edit the files in {args.export}/{top} as follows:')
for m in macros:
    print(f'- Remove USE_{m} line(s) in Makefile.in.')
if 'GSL' not in macros:
    print('- Uncomment lines in Makefile.in related to USE_GSL.')
if len(macros) > 0:
    print('- Also consider removing references to deleted macros in Makefile.')
if 'DEM' in macros:
    print('- Remove SSDEM-related files from other parts of the codebase.')
if 'SPRINGS' in macros:
    print('- Remove springs.h from include and src/pkdgrav directories.')
if 'WALLS' in macros:
    print('- Remove wallsio.c and wallsio.h from src/pkdgrav and src/ss.')
if 'ORIGIN_HISTOGRAM' in extras:
    print('- Remove ORIGIN_HISTOGRAM_BINS lines in Makefile.in.')
print('- It may also be necessary to alter etc/ss.par.')
print('- In addition, README and README.COPYRIGHT should be replaced.')
if 'COLLISIONS' in macros:
    print('- Removing COLLISIONS effectively means removing most of the\n'
          '  ss-related components from the code base. This will require\n'
          '  additional manual steps -- you are on your own!')
print(f'\nTo make a tarball, do the following in the {args.export} directory:')
print(f'   tar cvfz tarball-name.tar.gz {top}')
print(f'(Consider renaming "{top}" to include a version label\n'
      'or date and use the same scheme for the tarball name.)')
