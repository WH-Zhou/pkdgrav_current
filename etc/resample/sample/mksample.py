#!/usr/bin/env python3

import glob
import os
import subprocess

# specify ssgen arguments here...(rerun demparams if needed!)...

args = ["-r", "100.0",  # minimum particle radius in cm
        "-d", "2.0",  # particle mass density in g/cc
        "-p", "-3.5,3.0,2",  # power-law index, size ratio, distribution
        "-x", "25.0",  # x dimension of region in m
        "-y", "25.0",  # y dimension of region in m
        "-z", "25.0",  # z dimension of region in m
        "-e",  # to use ellipsoidal region
        "100"]  # number of particles to generate

print(args)
try:
    subprocess.run(["ssgen"] + args, check=True)
except subprocess.CalledProcessError:
    print("FATAL: Error running ssgen.")
    exit(1)
try:
    subprocess.run(["./pkdgrav", "ss.par"], check=True)
except subprocess.CalledProcessError:
    print("FATAL: Error running pkdgrav.")
    exit(1)
ssfiles = sorted(glob.glob("ss.[0-9]*[0-9]"))
os.symlink("ssgen.ss", "ss." + "0" * (len(ssfiles[0].split()[-1]) - 3))
try:
    subprocess.run("mkmov.py", check=True)
except subprocess.CalledProcessError:
    print("WARNING: Error running mkmov.py.")
print(f"Now use {ssfiles[-1]} as a replacement ss file for ssresample.py.")
