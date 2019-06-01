#!/usr/bin/python
import os

parameters = \
    """#!/bin/bash
#SBATCH --workdir /scratch/milenkov/run22/
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16000
#SBATCH --time 00:05:00

/home/milenkov/deneb/bsfft/optimized/precompute_filters -n %d -b %d -p /scratch/milenkov/filters10/
"""

folderName = "filter_sbatch_files"
if not os.path.exists(folderName):
    os.makedirs(folderName)

for n in range(10, 23):
    b = 2
    while b <= pow(2, n) / 4:
        file = open(folderName + ("/create_filter_%d_%d.run" % (n, b)), "w")
        file.write(parameters % (n, b))
        file.close()
        b *= 2
