import sys
import os
import numpy as np
import pyvista as pv
import ssio

# Constants
G = 6.67430e-11  # m^3 kg^-1 s^-2
yr = 3.15576e7   # seconds
au = 1.496e11    # m
M_sun = 1.989e30 # kg

def generate_filenames(path, start, end, step):
    return [f"{path}/ss.{str(i).zfill(10)}" for i in range(start, end, step)]

def ss_to_vtp(ss_file, output_file):
    data = ssio.read_SS(ss_file)[0]

    # Unpack and convert fields
    particle_id = data[0].astype(np.int32)
    index       = data[1].astype(np.int32)
    mass        = data[2] * M_sun
    radius      = data[3] * au
    x, y, z     = data[4]*au, data[5]*au, data[6]*au
    vx, vy, vz  = data[7]*au/yr*2*np.pi, data[8]*au/yr*2*np.pi, data[9]*au/yr*2*np.pi
    sx, sy, sz  = data[10]/yr*2*np.pi, data[11]/yr*2*np.pi, data[12]/yr*2*np.pi
    color       = data[13].astype(np.int32) if data.shape[0] > 13 else np.zeros_like(particle_id)

    # Create point cloud
    points = np.column_stack((
        np.asarray(x, dtype=np.float64),
        np.asarray(y, dtype=np.float64),
        np.asarray(z, dtype=np.float64)
    ))

    cloud = pv.PolyData(points)

    # Add scalar and vector fields
    cloud['particle_id'] = particle_id
    cloud['index']       = index
    cloud['mass']        = mass
    cloud['radius']      = radius
    cloud['color']       = color
    cloud['velocity']    = np.column_stack((vx, vy, vz)).astype(np.float64)
    cloud['spin']        = np.column_stack((sx, sy, sz)).astype(np.float64)

    cloud.save(output_file)

def write_pvd(path, start, end, step):
    pvd_path = os.path.join(path, "frames.pvd")
    with open(pvd_path, 'w') as f:
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <Collection>\n')

        for i in range((end - start) // step):
            frame_id = start + i * step
            filename = f"frame_{frame_id}.vtp"
            f.write(f'    <DataSet timestep="{frame_id}" group="" part="0" file="{filename}"/>\n')

        f.write('  </Collection>\n')
        f.write('</VTKFile>\n')

    print(f"✅ Wrote time series PVD: {pvd_path}")

def main():
    if len(sys.argv) != 5:
        print("Usage: python ss_to_vtp_batch_with_pvd.py <path> <start> <end> <step>")
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

    write_pvd(path, start, end, step)

if __name__ == "__main__":
    main()
