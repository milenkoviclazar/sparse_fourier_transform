#!/usr/bin/python
import os

parameters = \
    """#!/bin/bash
#SBATCH --workdir /scratch/milenkov/run_signals/
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16000
#SBATCH --time 00:01:30

/home/milenkov/deneb/bsfft/optimized/generate_signals -n %d -0 %d -1 %d -c 20 -p /scratch/milenkov/random_signals/
"""
# TODO: check --workdir
# TODO: check --time limit
# TODO: check cpp executable name
# TODO: check -p switch for the executable
# TODO: check folderName for '*.run' files
# TODO: ranges for the parameters
# TODO: make!!!
# TODO: generate the '*.run' file and try runnning the cpp executable with the generated parameters

folderName = "signal_sbatch_files"
if not os.path.exists(folderName):
    os.makedirs(folderName)


for n in range(10, 23):
    for k0 in range(1, 33):
        k1 = 4
        while k1 <= 128:
            file = open(folderName + ("/generate_signal_%d_%d_%d.run" % (n, k0, k1)), "w")
            file.write(parameters % (n, k0, k1))
            file.close()
            k1 *= 2
