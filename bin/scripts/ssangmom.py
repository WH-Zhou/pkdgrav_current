#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 09:48:16 2020

@author: dcr

Utility script to plot angular momentum evolution from ss files.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    import ssio
except ModuleNotFoundError:
    print("Error: ssio library not found.")
    print("Add location of ssio.py to PYTHONPATH environment variable.")
    print("Or copy ssio.py to this directory.")
    print("Look for ssio.py in the source code bin/scripts directory.")
    sys.exit(1)

parser = argparse.ArgumentParser(
    description='Plot angular momentum evolution from ss files.')
parser.add_argument('ssfiles', metavar='ssfile', nargs='+',
                    help='ss files to process.')
args = parser.parse_args()

# basic units...

units = {
    'mass': 1.989e30,  # 1 solar mass in kg (multiply to convert to kg)
    'length': 1.495e11,  # 1 au in m (convert to m)
    'time': 365.25 * 24 * 3600 / (2 * np.pi)  # 1 yr/2 pi in s (convert to s)
}

# derived units...

units['speed'] = units['length'] / units['time']  # convert to m/s
units['spin'] = 2 * np.pi / units['time']  # convert to rad/s
units['density'] = units['mass'] / units['length']**3  # convert to kg/m^3
units['ang mom'] = units['mass'] * units['speed'] * units['length']
units['energy'] = units['mass'] * units['speed'] / units['time']

units['1/h'] = 3600 * units['spin'] / (2 * np.pi)**2  # to 1/h

# read the data...

t, L = [], []
for ssfile in args.ssfiles:
    print(ssfile)
    try:
        header, data = ssio.read_SS(ssfile, 'y')
    except ssio.ReadError as msg:
        print('Read error:', msg)
        sys.exit(1)
    t.append(header[0])
    # handy shortcuts...
    m = data[2, :]
    R = data[3, :]
    pos = data[4:7, :]
    vel = data[7:10, :]
    spin = data[10:13, :]
    m_tot = sum(m)
    mom = 0.4*m*R**2
    assert m_tot > 0
    # first get center-of-mass data...
    com_pos = (m[np.newaxis, :] * pos).sum(axis=1) / m_tot
    com_vel = (m[np.newaxis, :] * vel).sum(axis=1) / m_tot
    print('Center of mass position [m]', com_pos * units['length'])
    print('Center of mass velocity [m/s]', com_vel * units['speed'])
    # subtract c-o-m pos & vel from each particle...
    pos -= com_pos[:, np.newaxis]
    vel -= com_vel[:, np.newaxis]
    # calculate total angular momentum...
    ang_mom = m[np.newaxis, :]*np.cross(pos, vel, axisa=0, axisb=0).T \
        + mom[np.newaxis, :]*spin
    ang_mom_tot = ang_mom.sum(axis=1)
    print('Total angular momentum [kg m^2/s]:', ang_mom_tot * units['ang mom'])
    L.append(np.linalg.norm(ang_mom.sum(axis=1)))

plt.plot(np.array(t)*units['time']/(24*3600), np.array(L)*units['ang mom'],
         '.-')
plt.xlabel('Time [d]')
plt.ylabel(r'Total Angular momentum [kg m$^2$ s$^{-1}$]')
plt.show()
