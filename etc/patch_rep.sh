#!/bin/bash

NPROC=16

echo "" > PATCH_CHECK.txt

if [ ! -e patchic.ss ]; then
    mpirun -n $NPROC ./patchic
    mpirun -n $NPROC ./pkdgrav +rejects patch.par | tee -a PATCH_CHECK.txt
fi

set -o pipefail

while :; do
    echo "@@ running patchic @@"
    mpirun -n $NPROC ./patchic -r 2>&1 | tee -a PATCH_CHECK.txt

    if [ $? -ne 0 ]; then
        echo $?
        echo "error in patchic!"
        exit
    fi

    echo "@@ running pkdgrav @@"
    mpirun -n $NPROC ./pkdgrav +rejects patch.par 2>&1 | tee -a PATCH_CHECK.txt

    if [ $? -ge 2 ]; then
        echo $?
        echo "error in pkdgrav!"
        exit
    fi

    echo "@@ checking @@"
    found=`grep "rejects\? found!" PATCH_CHECK.txt | tail -n 1`
    echo "$found" | tee -a PATCH_CHECK.txt
    if [ ! -z "$found" ]; then
        echo $found
        continue
    fi
    exit
done

