import os
import subprocess
import sys
from pathlib import Path

def generate_filenames(path, start, end, step):
    return [path + f"/ss.{str(i).zfill(10)}" for i in range(start, end, step)]

def run_draw_image(file):
    # print(f"🎨 Rendering image from: {file}")
    subprocess.run(["python3", "draw_image.py", file], check=True)

def create_file_list(image_files, duration, list_filename="file_list.txt"):
    # print("📝 Creating file_list.txt for ffmpeg...")
    with open(list_filename, 'w') as f:
        for img in image_files:
            f.write(f"file '{img}'\n")
            f.write(f"duration {duration}\n")
        # Repeat last file to avoid ffmpeg cutoff
        f.write(f"file '{image_files[-1]}'\n")

def create_video_from_list(output, list_file="file_list.txt"):
    # print("🎬 Generating video with ffmpeg...")
    subprocess.run([
        "ffmpeg", "-f", "concat", "-safe", "0",
        "-i", list_file,
        "-vsync", "vfr",
        "-pix_fmt", "yuv420p",
        "-c:v", "libx264",
        output
    ], check=True)
    # print(f"✅ Video created: {output}")

def cleanup_images(image_files):
    for img in image_files:
        try:
            os.remove(img)
        except OSError as e:
            print(f"Error deleting {img}: {e}")

def main(path, start_index, end_index, step=10000, duration=0.5, output = "output"):

    filenames = generate_filenames(path, start_index, end_index, step)

    rendered_images = []
    pov_images = []
    bases = []


    for file in filenames:
        run_draw_image(file)
        # base = Path(file).stem  # e.g., ss.0000001000
        base = os.path.basename(file)        # 'initcond.ss'
        if file.endswith(".ss") or file.endswith(".r"):
            base = os.path.splitext(base)[0]
        bases.append(base)
        rendered_images.append(f"{base}.png")
        pov_images.append(f"{base}.pov")

        print(rendered_images[-1], pov_images[-1])

    create_file_list(rendered_images, duration)
    create_video_from_list(output=output)
    cleanup_images(rendered_images)
    cleanup_images(pov_images)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit(1)

    path = str(sys.argv[1])
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    step = int(sys.argv[4])
    duration = float(sys.argv[5])
    output = str(sys.argv[6] + ".mp4")
    main(path, start, end, step, duration, output)
