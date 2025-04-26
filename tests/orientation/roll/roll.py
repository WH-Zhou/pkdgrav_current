import glob
from ssio import *
import numpy as np
import matplotlib.pyplot as plt

figsize = (8, 8)

files = sorted(glob.glob("ss.[0-9]*[0-9]"))

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

t, r, v, w = np.array([]), np.array([]), np.array([]), np.array([])

for file in files:
    header, data = read_SS(file, "yes")
    assert header[1] == 1
    t = np.append(t, header[0])
    r = np.append(r, [data[4, 0], data[5, 0], data[6, 0]])
    v = np.append(v, [data[7, 0], data[8, 0], data[9, 0]])
    w = np.append(w, [data[10, 0], data[11, 0], data[12, 0]])

t = t * units['time']
r = r * units['length'] * 100  # cm
v = v * units['speed'] * 100  # cm/s
w = w * units['1/h']

r.shape = (len(t), 3)
v.shape = (len(t), 3)
w.shape = (len(t), 3)

plt.figure(figsize=figsize)
plt.suptitle("Particle Rolling in $+x$-Direction")
plt.subplot(311)
plt.plot(t, r[:, 0], '.-')
plt.xticks(visible=False)
plt.ylabel("$x$ [cm]")
plt.grid()
plt.subplot(312)
plt.plot(t, v[:, 0], '.-')
plt.xticks(visible=False)
plt.ylabel("$v_x$ [cm/s]")
plt.grid()
plt.subplot(313)
plt.plot(t, w[:, 1], '.-')
plt.xlabel("$t$ [s]")
plt.ylabel(r"$\omega_y$ [h$^{-1}$]")
plt.grid()
plt.subplots_adjust(hspace=0)
plt.show()
