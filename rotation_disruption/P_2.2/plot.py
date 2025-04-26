import pkdtools as pt
import numpy as np
import matplotlib.pyplot as plt

def plot_file(filename):
    data = pt.ss_in(filename, units='mks')
    x = np.zeros(len(data))
    y = np.zeros(len(data))
    z = np.zeros(len(data))
    for i in range(len(data)):
        x[i] = data[i].x
        y[i] = data[i].y
        z[i] = data[i].z
    plt.plot(x,y,'.')
    plt.xlim(-10e4,10e4)
    plt.ylim(-10e4,10e4)
    plt.pause(0.1)
    plt.cla()
    # plt.show()

def make_animation():
    max_index = 2000000
    min_index = 10000
    delta_index = 10000
    for i in range(min_index, max_index, delta_index):
        filename = 'ss.'+ '0'*(7-len(str(i)))  +  str(i)
        # print(filename)
        plot_file(filename)
    # plot_file('ss.0')

fig = plt.figure()
make_animation()