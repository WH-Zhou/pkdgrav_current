import glob
import numpy as np
import matplotlib.pyplot as plt
from ssio import *

t_conv = 365.25 * 24 * 3600 / (2 * np.pi)  # pkdgrav spin units to seconds
w_conv = 1 / t_conv  # pkdgrav time units to seconds

figsize = (8, 8)

files = sorted(glob.glob("ss.[0-9]*[0-9]"))

t = np.array([])
w = np.array([])

for file in files:
    print(file)
    header, data = read_SS(file, "yes")
    assert header[1] == 2  # two particles
    t = np.append(t, header[0] * t_conv)
    w = np.append(w,
                  np.array([data[10, 0], data[11, 0], data[12, 0],
                            data[10, 1], data[11, 1], data[12, 1]]) * w_conv)

w.shape = (len(t), 2, 3)

x0 = w[:, 0, 0]
y0 = w[:, 0, 1]
z0 = w[:, 0, 2]
x1 = w[:, 1, 0]
y1 = w[:, 1, 1]
z1 = w[:, 1, 2]

# plt.plot(t, z0, t, z1)

plt.figure(figsize=figsize)
plt.suptitle("Particle 0")
plt.plot(t, x0, '.-', label="x")
plt.plot(t, y0, '.-', label="y")
plt.plot(t, z0, '.-', label="z")
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Spin (rad/s)")
plt.figure(figsize=figsize)
plt.suptitle("Particle 1")
plt.plot(t, x1, '.-', label="x")
plt.plot(t, y1, '.-', label="y")
plt.plot(t, z1, '.-', label="z")
plt.legend()
plt.xlabel("Time (s)")
plt.ylabel("Spin (rad/s)")

# mag0 = [np.linalg.norm(x) for x in w[:, 0]]
# mag1 = [np.linalg.norm(x) for x in w[:, 1]]
# plt.plot(t, mag0, 'x-', label="0")
# plt.plot(t, mag1, 'o-', label="1")
# plt.legend()

plt.show()
