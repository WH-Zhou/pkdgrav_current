import os
import subprocess
import sys
from pathlib import Path

def generate_filenames(start, end, step):
    return [f"ss.{str(i).zfill(10)}" for i in range(start, end, step)]

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

def main(start_index, end_index, step=10000, duration=0.5, output = "output"):
    filenames = generate_filenames(start_index, end_index, step)

    rendered_images = []
    pov_images = []
    for file in filenames:
        run_draw_image(file)
        # base = Path(file).stem  # e.g., ss.0000001000
        rendered_images.append(f"{file}.png")
        pov_images.append(f"{file}.pov")

    create_file_list(rendered_images, duration)
    create_video_from_list(output=output)
    cleanup_images(rendered_images)
    cleanup_images(pov_images)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(1)

    start = int(sys.argv[1])
    end = int(sys.argv[2])
    step = int(sys.argv[3])
    duration = float(sys.argv[4])
    output = str(sys.argv[5] + ".mp4")
    main(start, end, step, duration, output)
