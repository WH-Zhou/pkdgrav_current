#!/usr/bin/env python3

import argparse
import random
import numpy as np

try:
    import ssio
except ModuleNotFoundError:
    print('Error: ssio library not found.')
    print('Add location of ssio.py to PYTHONPATH environment variable.')
    print('Or copy ssio.py to this directory.')
    print('Look for ssio.py in the source code bin/scripts directory.')
    exit(1)


# helper function for reading ss data (via ssio.py)...

def read_data(ssfile, header=True, recenter=True):
    # header format:
    #  header[0] = time
    #  header[1] = number of particles
    #  header[2] = magic number (ignore)
    # data format (per particle):
    #  data[0,:]     = order (ignore)
    #  data[1,:]     = original index
    #  data[2,:]     = mass
    #  data[3,:]     = radius
    #  data[4:7,:]   = position (3-vector)
    #  data[7:10,:]  = velocity (3-vector)
    #  data[10:13,:] = spin (3-vector)
    #  data[13,:]    = color (integer)
    try:
        if header:
            head, data = ssio.read_SS(ssfile, 'yes')
        else:
            data = ssio.read_SS(ssfile)
    except ssio.ReadError as msg:
        print(f'Error reading from "{ssfile}":', msg)
        exit(1)
    if recenter:
        total_mass = np.sum(data[2, :])
        assert total_mass > 0.0
        com = np.sum(data[2, :] * data[4:7, :], axis=1) / total_mass
        data[4:7, :] -= com[:, np.newaxis]  # broadcast vector across array
        com = np.sum(data[2, :] * data[7:10, :], axis=1) / total_mass
        data[7:10, :] -= com[:, np.newaxis]
    if header:
        return head, data
    else:
        return data


# handle command-line arguments...

parser = argparse.ArgumentParser(
        description='Replace particles in ss file with groups of particles.')
parser.add_argument('ssfile', type=str, help='source ssfile')
parser.add_argument('-k', '--keep', action='store_true',
                    help='keep particles that are not replaced')
parser.add_argument('-o', '--outfile', type=str, default='resampled.ss',
                    help='name of output file (default %(default)s)')
parser.add_argument('-p', '--particles', type=str, default='particles.txt',
                    help='text file listing particles to replace (or "all"; '
                    'default %(default)s)')
parser.add_argument('-r', '--replacements', type=str,
                    default='replacements.txt',
                    help='text file listing filenames of replacement ss files '
                    '(default %(default)s)')
parser.add_argument('-s', '--sepmult', type=float, default=1.01,
                    help='separation multiplier for replacement particles '
                    '(default %(default)g)')
args = parser.parse_args()
if args.sepmult <= 0.0:
    print('Separation multiplier must be positive.')
    exit(1)

# read source ss file...

src_head, src_data = read_data(args.ssfile, recenter=False)

# set/read particle list...

if args.particles == 'all':
    particles = range(src_head[1])
else:
    try:
        particles = sorted(
                [int(x) for x in open(args.particles, 'r').readlines()])
    except IOError:
        print(f'Unable to access particle list in "{args.particles}".')
        exit(1)

# get list of replacements...

try:
    replacements = \
        [x.strip() for x in open(args.replacements, 'r').readlines()]
except IOError:
    print(f'Unable to read list of replacements in "{args.replacements}".')
    exit(1)

# carry out replacements...

new_data = np.empty((14, 0))

old_idx = 0
for particle in particles:
    # concatenate source data to output, up to next replacement particle...
    if args.keep:
        new_data = np.c_[new_data, src_data[:, old_idx:particle]]
    old_idx = particle
    random.shuffle(replacements)
    try:
        rep_data = read_data(replacements[0], header=False)
    except IOError:
        print(f'Unable to access ss data in "{replacements[0]}".')
        exit(1)
    # get ratio of original particle radius to sphere enclosing replacements...
    scale = src_data[3, old_idx] / \
        np.max([np.linalg.norm(r) for r in rep_data[4:7, :].T])
    # set properties of replacement particles...
    rep_data[1, :] = src_data[1, old_idx]  # inherit original index
    rep_data[2, :] *= scale**3  # mass goes as length scale cubed
    rep_data[3, :] *= scale  # radius goes as length scale
    new_pos = rep_data[4:7, :] * scale * args.sepmult
    rep_data[4:7, :] = src_data[4:7, old_idx][:, np.newaxis] + new_pos
    wxr = np.array([np.cross(src_data[10:13, old_idx], r) for r in new_pos.T])
    rep_data[7:10, :] += src_data[7:10, old_idx][:, np.newaxis] + wxr.T
    rep_data[10:13, :] = src_data[10:13, old_idx][:, np.newaxis]
    rep_data[13, :] = src_data[13, old_idx]  # inherit color
    new_data = np.c_[new_data, rep_data]
    old_idx += 1

# concatenate remaining source data to output, if desired...
if args.keep:
    new_data = np.c_[new_data, src_data[:, old_idx:-1]]

try:
    ssio.write_SS(new_data, args.outfile, time=src_head[0])
except ssio.WriteError as msg:
    print(f'Error writing to "{args.outfile}":', msg)
