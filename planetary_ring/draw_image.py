import subprocess
import os
import sys

def run_command(command):
    """Utility to run a shell command and print its output."""
    # print(f"Running: {' '.join(command)}")
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error:\n{result.stderr}")
        sys.exit(result.returncode)
    print(result.stdout)

def generate_image(filename):
    # Step 1: Run ssdraw
    path = 'results/'
    run_command(["ssdraw", path + filename])

    # Assuming output is filename.pov.gz
    base = filename
    if filename.endswith(".ss") or filename.endswith(".r"):
        base = os.path.splitext(filename)[0]
    gz_file = f"{base}.pov.gz"
    pov_file = f"{base}.pov"
    output_image = f"{base}.png"


    # Step 2: Decompress with gunzip
    run_command(["gunzip", gz_file])

    # Step 3: Render with POV-Ray
    run_command([
        "povray",
        f"+I{pov_file}",
        f"+O{output_image}",
        "+W800",
        "+H600",
        "+FN",     # PNG output
        "+A",      # Anti-aliasing
        "+Q11"     # Quality level
    ])

    # print(f"✅ Image generated: {output_image}")

# --- Usage ---
if __name__ == "__main__":
    if len(sys.argv) != 2:
        # print("Usage: python render_pipeline.py <filename>")
        sys.exit(1)
    input_filename = sys.argv[1]
    generate_image(input_filename)
