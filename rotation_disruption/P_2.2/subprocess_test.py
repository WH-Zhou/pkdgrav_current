import sys
sys.path.append('/scratch/wenzhou/pkdgrav_current_master/bin/scripts')

import subprocess
import pkdtools as pt
import numpy as np

def reset_initcond(init_period, file1, file2):
    # init_period: h
    data = pt.ss_in(file1, units='mks')
    period = data.period()     # s
    # print(period)
    ratio = period /(init_period * 3600)
    for x in data:
        x.vx = x.vx * ratio
        x.vy = x.vy * ratio
        x.vz = x.vz * ratio
    new_data = pt.Assembly(*[x for x in data], units='mks')
    pt.ss_out(new_data, file2)
    # print(new_data.period()/3600)


# def plot_file(filename):
#     data = pt.ss_in(filename, units='mks')
#     x = np.zeros(len(data))
#     y = np.zeros(len(data))
#     z = np.zeros(len(data))
#     for i in range(len(data)):
#         x[i] = data[i].x
#         y[i] = data[i].y
#         z[i] = data[i].z
#     plt.plot(x,y,'.')
#     plt.show()

def get_displacement(init, final):
    init_data = pt.ss_in(init, units='mks')
    final_data = pt.ss_in(final, units='mks')
    displacement = np.zeros([len(init_data),4])
    for i in range(len(init_data)):
        displacement[i,0] = final_data[i].x - init_data[i].x
        displacement[i,1] = final_data[i].y - init_data[i].y
        displacement[i,2] = final_data[i].z - init_data[i].z
        displacement[i,3] = np.sqrt(displacement[i,0]**2 + displacement[i,1]**2 + displacement[i,2]**2)
    return displacement


# Open the subprocess and interact with it


while True:
    p1 = subprocess.run(['sbatch', '/scratch/wenzhou/pkdgrav_current_master/rotation_disruption/P_2.2/runJob.sh'])
    p2 = subprocess.run(['sleep', '43200'])
    displacement = get_displacement('initcond.ss', 'ss.2000000')
    print('Maximum displacement')
    if max(displacement[:,3]) > 1:
        reset_initcond(2.2, 'ss.2000000', 'initcond.ss')
    else:
        break

print('done')


# Send inputs and read outputs as necessary
# proc.stdin.write("\n".encode())  # Simulate pressing ENTER
# proc.stdin.write("\n".encode())

# proc.stdin.write("ss.2000000\n".encode())
# proc.stdin.write("7\n".encode()) # Send "3" and then simulate ENTER
# proc.stdin.write("\n".encode())
# proc.stdin.write("\n".encode())
# proc.stdin.write("\n".encode())
# proc.stdin.write("0 0 0.5\n".encode())
# proc.stdin.write("0\n".encode())
# proc.stdin.write("\n".encode())
# proc.stdin.write("\n".encode())
# proc.stdin.write("initcond.ss\n".encode())
# proc.stdin.write("y\n".encode())
# Ensure you close the stdin to notify the process that you're done sending data

# proc.stdin.write('\n\nss.2000000\n7\n\n\n\n0 0 0.5\n0\n\n\ninitcond.ss\ny\n'.encode())


# proc.stdin.close()

# print('done')

# Wait for the process to complete and get the output
#output, errors = proc.communicate()

# Process the output and errors as needed
#print(output)

# Check the return code to ensure the process ran successfull#if proc.returncode != 0:
#    print(f"Error: {errors}")

