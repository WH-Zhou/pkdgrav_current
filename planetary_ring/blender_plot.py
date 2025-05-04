import bpy
import csv

csv_path = "/Users/wenzhou/Documents/pkdgrav_current_master/planetary_ring/initcond.csv" # Update this

# Clear objects
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

with open(csv_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for i, row in enumerate(reader):
        x, y, z = float(row['x']), float(row['y']), float(row['z'])
        r = float(row['r'])
        bpy.ops.mesh.primitive_uv_sphere_add(radius=r, location=(x, y, z))
        bpy.context.active_object.name = f"Particle_{i}"


