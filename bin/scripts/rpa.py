#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

# global customizations...enable as desired...
# import matplotlib as mpl
# mpl.rcParams['axes.formatter.limits'] = (-4, 4)
# mpl.rcParams['axes.formatter.useoffset'] = False

figsize = (12, 8)  # 12'x8' (may depend on the display!)

parser = argparse.ArgumentParser(description='Analyze and plot rpa data.')
parser.add_argument('datafile', default='rpa.out', nargs='?',
                    help='name of rpa data file (default: "%(default)s")')
parser.add_argument('--nodisplay', help='disable output of plots to screen',
                    action='store_true')
parser.add_argument('--nofiles', help='disable output of png files',
                    action='store_true')
parser.add_argument('--nopause', help='disable pausing between plots',
                    action='store_true')
parser.add_argument('--nostats', help='disable stats output',
                    action='store_true')
group = parser.add_mutually_exclusive_group()
group.add_argument('--dotsonly', help='only plot points, not lines',
                   action='store_true', default=False)
group.add_argument('--linesonly', help='only plot lines, not points',
                   action='store_true', default=False)
group.add_argument('--dotsandlines',
                   help='plot both points and lines (default)',
                   action='store_true', default=True)
args = parser.parse_args()
if args.nodisplay:
    args.nopause = True
if args.nodisplay and args.nofiles and args.nostats:
    print('Nothing to do!')
    exit(0)
if args.dotsonly or args.linesonly:
    args.dotsandlines = False

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
units['ang mom'] = units['mass'] * units['length'] * units['speed']
units['energy'] = units['mass'] * units['speed']**2

units['1/h'] = 3600 * units['spin'] / (2 * np.pi)**2  # to 1/h

# read the data...

try:
    data = np.genfromtxt(args.datafile, dtype=None, names=True, skip_header=2,
                         encoding='UTF-8')
    # handle special case of just 1 time value (force the data to an array)...
    if len(data.shape) == 0:
        data = np.array([data])
    if len(data) == 0:
        print('No particle data found in rpa.out (check the file).')
        exit(1)
except IOError:
    print(f'Unable to open "{args.datafile}"...use "-h" for help.')
    exit(1)

# convenient aliases for data access... (column labels in file header)
# NOTE: for backwards compatability, if new data columns are added to
# this script in the future, use try/except to skip over missing data
# in older formats.
time = data[:]['time'] * units['time'] / (24 * 3600)  # days
world_radius = data[:]['world_radius']
world_pos = data[:]['world_pos']
world_vel = data[:]['world_vel']
world_spin = data[:]['world_spin']
max_mass = data[:]['max_mass']
max_vel = data[:]['max_vel']
max_spinx = data[:]['max_spinx']
max_spiny = data[:]['max_spiny']
max_spinz = data[:]['max_spinz']
eff_spin = data[:]['eff_spin']
rot_idx = data[:]['rot_idx']
ang_mom = data[:]['ang_mom']
ke = data[:]['ke']
max_major = data[:]['max_major']
max_inter = data[:]['max_inter']
max_minor = data[:]['max_minor']
max_den = data[:]['max_den']
n1 = data[:]['n1']
n2 = data[:]['n2']
n = data[:]['n']
m1 = data[:]['m1']
m2 = data[:]['m2']
m = data[:]['m']
vd1 = data[:]['vd1']
vd2 = data[:]['vd2']
vd = data[:]['vd']
mAcc = data[:]['mAcc']
mOrb = data[:]['mOrb']
mEsc = data[:]['mEsc']

# derived quantities...

mTot = m1 + m2 + m  # total mass
try:
    q2 = max_inter / max_major  # intermediate axis ratio
except ValueError:
    print('Zero length major axis?')
    exit(1)
q3 = max_minor / max_major  # minor axis ratio
e = 1 - 0.5 * (q2 + q3)  # ellipticity


def pause(no_plots_yet=False):
    if args.nopause:
        return
    if not no_plots_yet and not args.nodisplay:
        plt.show(block=False)
    input('Press ENTER to continue...')
    if not no_plots_yet:
        plt.close()


if not args.nostats:

    # print a bunch of stats...

    def display(label, values, unit):
        print('%s = %g %s%s(avg %g +/- %g)' %
              (label, values[-1], unit, '' if unit == '' else ' ',
               np.mean(values), np.std(values)))

    print('Final stats (at time t = %g d):' % time[-1])
    try:
        display('Maximum expansion', world_radius / world_radius.min(),
                'times')
    except ValueError:
        print('Non-positive maximum expansion...no particle separation?')
        exit(1)
    display('Bary. distance', world_pos * units['length'] / 1000, 'km')
    display('Bary. speed', world_vel * units['speed'], 'm/s')
    display('Bulk spin', world_spin * units['1/h'], '1/h')
    print('Final stats of biggest rubble pile:')
    try:
        display('Mass', max_mass / mTot, '* total')
    except ValueError:
        print('Non-positive total mass...no particle mass?')
        exit(1)
    display('X-spin', max_spinx * units['1/h'], '1/h')
    display('Y-spin', max_spiny * units['1/h'], '1/h')
    display('Z-spin', max_spinz * units['1/h'], '1/h')
    display('Effective spin', eff_spin * units['1/h'], '1/h')
    display('Rotation index', rot_idx, '')
    display('Ang. mom. mag.', ang_mom * units['ang mom'], 'kg m^2/s')
    display('Kinetic energy', ke * units['energy'], 'J')
    display('Major axis', max_major * units['length'] / 1000, 'km')
    display('Intermediate axis', max_inter * units['length'] / 1000, 'km')
    display('Minor axis', max_minor * units['length'] / 1000, 'km')
    display('Density', max_den * units['density'] / 1000, 'g/cc')
    print('Final grouping stats:')
    display('1-particle', n1, 'group%s' % ('' if n1[-1] == 1 else 's'))
    display('2-particle', n2, 'group%s' % ('' if n2[-1] == 1 else 's'))
    display('Rubble piles', n, 'group%s' % ('' if n[-1] == 1 else 's'))
    display('1-particle mass', m1 / mTot, '* total')
    display('2-particle mass', m2 / mTot, '* total')
    display('Rubble piles mass', m / mTot, '* total')
    display('1-particle vel. disp.', vd1 * units['speed'], 'm/s')
    display('2-particle vel. disp.', vd2 * units['speed'], 'm/s')
    display('Rubble piles vel. disp.', vd * units['speed'], 'm/s')
    display('Accreting mass', mAcc / mTot, '* total')
    display('Orbiting mass', mOrb / mTot, '* total')
    display('Escaping mass', mEsc / mTot, '* total')
    pause(no_plots_yet=True)


def style(color, symbol, line):
    style_str = color
    if args.dotsonly or args.dotsandlines:
        style_str += symbol
    if args.linesonly or args.dotsandlines:
        style_str += line
    return style_str


# plot 'world' radius vs. time...

plt.figure(figsize=figsize)
plt.plot(time, world_radius / world_radius.min(), style('m', 'x', '-'))
plt.xlabel('Time (d)')
plt.ylabel('Maximum Expansion')
plt.grid()
plt.suptitle('Maximum Relative Distance Between Particles')
if not args.nofiles:
    plt.savefig('rpa_maxexp.png')
pause()

# world position, speed, spin vs. time...

plt.figure(figsize=figsize)
plt.subplot(311)
plt.plot(time, world_pos * units['length'] / 1000, style('r', 'x', '-'))  # km
plt.xticks(visible=False)
plt.ylabel('Barycenter Distance (km)')
plt.grid(axis='x')
plt.subplot(312)
plt.plot(time, world_vel * units['speed'], style('g', 'x', '-'))
plt.xticks(visible=False)
plt.ylabel('Barycenter Speed (m s$^{-1}$)')
plt.grid(axis='x')
plt.subplot(313)
plt.plot(time, world_spin * units['1/h'], style('b', 'x', '-'))  # 1/h
plt.xlabel('Time (d)')
plt.ylabel('Bulk Spin (h$^{-1}$)')
plt.grid(axis='x')
plt.subplots_adjust(hspace=0)
plt.suptitle('Bulk Distance, Speed, and Spin')
if not args.nofiles:
    plt.savefig('rpa_bulk.png')
pause()

# mass, speed, and bulk density of largest rubble pile...

plt.figure(figsize=figsize)
plt.subplot(311)
plt.plot(time, max_mass / mTot, style('r', 'x', '-'))
plt.xticks(visible=False)
plt.ylabel(r'Mass ($m_{\mathrm{tot}}$)')
plt.grid(axis='x')
plt.subplot(312)
plt.plot(time, max_vel * units['speed'], style('g', 'x', '-'))
plt.xticks(visible=False)
plt.ylabel('Speed (m s$^{-1}$)')
plt.grid(axis='x')
plt.subplot(313)
plt.plot(time, max_den * units['density'] / 1000, style('b', 'x', '-'))  # g/cc
plt.xlabel('Time (d)')
plt.ylabel('Bulk Density (g cm$^{-3}$)')
plt.grid(axis='x')
plt.subplots_adjust(hspace=0)
plt.suptitle('Mass, Speed, and Bulk Density of Largest Rubble Pile')
if not args.nofiles:
    plt.savefig('rpa_mvd.png')
pause()

# axis & spin evolution of largest rubble pile...

plt.figure(figsize=figsize)
plt.subplot(221)
plt.plot(time, max_major * units['length'] / 1000, style('r', 'x', '-'),
         label='$a_1$')  # km
plt.plot(time, max_inter * units['length'] / 1000, style('g', 'o', ':'),
         label='$a_2$')  # km
plt.plot(time, max_minor * units['length'] / 1000, style('b', '*', '--'),
         label='$a_3$')  # km
plt.xlabel('Time (d)')
plt.ylabel('Semiaxis Lengths $a_1$, $a_2$, $a_3$ (km)')
plt.legend()
plt.subplot(222)
plt.plot(time, e, style('r', 'x', '-'), label=r'$\epsilon$')
plt.plot(time, q2, style('g', 'o', ':'), label='$q_2$')
plt.plot(time, q3, style('b', '*', '--'), label='$q_3$')
plt.xlabel('Time (d)')
plt.ylabel(r'Ellipticity $\epsilon$ & Axis Ratios $q_2$, $q_3$')
plt.legend()
plt.subplot(223)
plt.plot(time, max_spinx * units['1/h'], style('r', 'x', '-'),
         label=r'$\omega_x$')
plt.plot(time, max_spiny * units['1/h'], style('g', 'o', ':'),
         label=r'$\omega_y$')
plt.plot(time, max_spinz * units['1/h'], style('b', '*', '--'),
         label=r'$\omega_z$')
plt.xlabel('Time (d)')
plt.ylabel(r'Spin Components $\omega_x$, $\omega_y$, $\omega_z$ (h$^{-1}$)')
plt.legend()
plt.subplot(224)
plt.plot(time, np.sqrt(max_spinx**2 + max_spiny**2 + max_spinz**2) *
         units['1/h'], style('m', 'x', '-'))
plt.xlabel('Time (d)')
plt.ylabel(r'Spin Magnitude $\omega$ (h$^{-1}$)')
plt.suptitle('Axis and Spin Evolution of Largest Rubble Pile')
if not args.nofiles:
    plt.savefig('rpa_aw.png')
pause()

# rotation properties of largest rubble pile...

plt.figure(figsize=figsize)
plt.subplot(221)
plt.plot(time, eff_spin * units['1/h'], style('m', 'x', '-'))
plt.xlabel('Time (d)')
plt.ylabel('Effective Spin (h$^{-1}$)')
plt.subplot(222)
plt.plot(time, rot_idx, style('m', 'x', '-'))
plt.xlabel('Time (d)')
plt.ylabel('Rotation Index')
plt.subplot(223)
plt.plot(time, ang_mom * units['ang mom'], style('m', 'x', '-'))
plt.xlabel('Time (d)')
plt.ylabel('Ang. Mom. Mag. (kg m$^2$ s$^{-1}$)')
plt.subplot(224)
plt.plot(time, ke * units['energy'], style('m', 'x', '-'))
plt.xlabel('Time (d)')
plt.ylabel('Kinetic Energy (J)')
plt.suptitle('Rotation Properties of Largest Rubble Pile')
if not args.nofiles:
    plt.savefig('rpa_rot.png')
pause()

# group stats...

plt.figure(figsize=figsize)
plt.subplot(311)
plt.plot(time, n1, style('r', 'x', '-'), label='$n_1$')
plt.plot(time, n2, style('g', 'o', ':'), label='$n_2$')
plt.plot(time, n, style('b', '*', '--'), label='$n$')
plt.xticks(visible=False)
plt.ylabel('No. Groups')
plt.grid(axis='x')
plt.legend()
plt.subplot(312)
plt.plot(time, m1 / mTot, style('r', 'x', '-'), label='$m_1$')
plt.plot(time, m2 / mTot, style('g', 'o', ':'), label='$m_2$')
plt.plot(time, m / mTot, style('b', '*', '--'), label='$m$')
plt.xticks(visible=False)
plt.ylabel(r'Mass ($m_{\mathrm{tot}}$)')
plt.grid(axis='x')
plt.legend()
plt.subplot(313)
plt.plot(time, vd1 * units['speed'], style('r', 'x', '-'), label=r'$\sigma_1$')
plt.plot(time, vd2 * units['speed'], style('g', 'o', ':'), label=r'$\sigma_2$')
plt.plot(time, vd * units['speed'], style('b', '*', '--'), label=r'$\sigma$')
plt.xlabel('Time (d)')
plt.ylabel('Vel. Disp. (m s$^{-1}$)')
plt.grid(axis='x')
plt.legend()
plt.subplots_adjust(hspace=0)
plt.suptitle('Group Stats')
if not args.nofiles:
    plt.savefig('rpa_gstats.png')
pause()

# mass groups...

plt.figure(figsize=figsize)
plt.plot(time, max_mass / mTot, style('m', 'x', '-'),
         label=r'$m_{\mathrm{max}}$')
plt.plot(time, mAcc / mTot, style('r', 'x', '-'),
         label=r'$m_{\mathrm{acc}}$')
plt.plot(time, mOrb / mTot, style('g', 'o', ':'),
         label=r'$m_{\mathrm{orb}}$')
plt.plot(time, mEsc / mTot, style('b', '*', '--'),
         label=r'$m_{\mathrm{esc}}$')
plt.xlabel('Time (d)')
plt.ylabel(r'$m_{\mathrm{max}}$, $m_{\mathrm{acc}}$, $m_{\mathrm{orb}}$, '
           r'$m_{\mathrm{esc}}$ ($m_{\mathrm{tot}}$)')
plt.grid()
plt.legend()
plt.suptitle('Mass Stats')
if not args.nofiles:
    plt.savefig('rpa_mstats.png')
pause()

if args.nopause and not args.nodisplay:
    plt.show()
