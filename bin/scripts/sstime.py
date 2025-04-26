#!/usr/bin/env python3
#
# Simple python script to set the time in an ss file header.  Uses ssio.py.
# DCR 4/8/16
#

import numpy as np
import argparse

try:
    import ssio
except ModuleNotFoundError:
    print('Error: ssio library not found.')
    print('Add location of ssio.py to PYTHONPATH environment variable.')
    print('Or copy ssio.py to this directory.')
    print('Look for ssio.py in the source code bin/scripts directory.')
    exit(1)

# converts system time units to seconds (multiply)...
sys_to_sec = 365.25 * 24 * 3600 / (2 * np.pi)

parser = argparse.ArgumentParser(description='Set time in ss file header.',
                                 epilog='Example: %(prog)s initcond.ss 0')
parser.add_argument('ssfile', help='name of ssfile')
parser.add_argument('time', nargs='?', type=float,
                    help='time in seconds or system units')
parser.add_argument('-q', '--quiet', action='store_true',
                    help='quiet output mode', dest='quiet')
parser.add_argument('-u', '--simunits', action='store_true',
                    help='use system units', dest='simunits')
args = parser.parse_args()

try:
    header, data = ssio.read_SS(args.ssfile, 'yes')
except ssio.ReadError as msg:
    print('Read error:', msg)
    exit(1)

if not args.quiet:
    t = header[0]
    print(f'{args.ssfile}: old time = {t} ({t * sys_to_sec} s).')

while args.time is None:
    try:
        units = 'system units' if args.simunits else 'seconds'
        args.time = float(input(f'Enter new time in {units}: '))
    except ValueError:
        print('Invalid time, try again.')

if not args.simunits:
    # convert from seconds to system time units...
    args.time = args.time / sys_to_sec

try:
    ssio.write_SS(data, args.ssfile, time=args.time)
except ssio.WriteError as msg:
    print('Write error:', msg)
    exit(1)

if not args.quiet:
    t = args.time
    print(f'{args.ssfile}: new time = {t} ({t * sys_to_sec} s).')
