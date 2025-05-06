import sys
import os
import numpy as np
import pyvista as pv
import ssio

G = 6.67430e-11  # m^3 kg^-1 s^-2
yr = 3.15576e7  # seconds
au = 1.496e11  # m
M_sun = 1.989e30  # kg


M_center = 6.39e23 # kg: Martian mass
r_center = 3.3895e6  # m: Martian radius
orbit_radius = 3 * r_center
omega_orb = np.sqrt(G * M_center / (orbit_radius)**3)

def generate_filenames(path, start, end, step):
    return [path + f"/ss.{str(i).zfill(10)}" for i in range(start, end, step)]

def ss_to_vtp(ss_file, output_file):
    # print(ss_file)
    data = ssio.read_SS(ss_file)[0]

    # Unpack fields
    particle_id = data[0].astype(np.int32)
    index       = data[1].astype(np.int32)
    mass        = data[2] * M_sun
    radius      = data[3] * au
    x, y, z     = data[4]* au, data[5]* au, data[6]* au
    vx, vy, vz  = data[7]* au/yr * 2 *np.pi, data[8]* au/yr * 2 *np.pi, data[9]* au/yr * 2 *np.pi
    sx, sy, sz  = data[10]/yr * 2 *np.pi, data[11]/yr * 2 *np.pi, data[12]/yr * 2 *np.pi + omega_orb
    color       = data[13].astype(np.int32) if data.shape[0] > 13 else np.zeros_like(particle_id)

    # Build point cloud
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    points = np.column_stack((x, y, z))  # shape (N, 3), dtype float64

    cloud = pv.PolyData(points)

    # Add fields
    cloud['particle_id'] = particle_id
    cloud['index']       = index
    cloud['mass']        = mass
    cloud['radius']      = radius
    
    vx = np.asarray(vx, dtype=np.float64)
    vy = np.asarray(vy, dtype=np.float64)
    vz = np.asarray(vz, dtype=np.float64)
    sx = np.asarray(sx, dtype=np.float64)
    sy = np.asarray(sy, dtype=np.float64)
    sz = np.asarray(sz, dtype=np.float64)

    cloud['velocity'] = np.column_stack((vx, vy, vz))
    cloud['spin'] = np.column_stack((sx, sy, sz))

    cloud['color']       = color

    # Save
    cloud.save(output_file)

def main():
    if len(sys.argv) != 5:
        print("Usage: python ss_to_vtp_batch.py <path> <start> <end> <step>")
        sys.exit(1)

    path  = str(sys.argv[1])
    start = int(sys.argv[2])
    end   = int(sys.argv[3])
    step  = int(sys.argv[4])

    ss_files = generate_filenames(path, start, end, step)

    for i, ss_file in enumerate(ss_files):
        frame_id = start + i * step
        output_file = os.path.join(path, f"frame_{frame_id}.vtp")
        print(f"Processing {ss_file} → {output_file}")
        ss_to_vtp(ss_file, output_file)

if __name__ == "__main__":
    main()
