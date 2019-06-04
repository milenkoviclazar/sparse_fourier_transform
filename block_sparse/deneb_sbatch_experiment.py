#!/usr/bin/python
import os
from datetime import datetime

timestamp_str = str(datetime.now()).replace(" ", "_")
scratchFolderName =  "/scratch/milenkov/" + timestamp_str + "/"

parameters = \
    """#!/bin/bash
#SBATCH --workdir %s
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 50GB
#SBATCH --time 06:00:00

/home/milenkov/deneb/sparse_fourier_transform/block_sparse/run_experiment -n %d -0 %d -1 %d -c 20 -x %f -y %f -z %f \
-p /scratch/milenkov/signals/ \
-q /scratch/milenkov/filters/
"""


folderName = "experiment_sbatch_files"
if not os.path.exists(folderName):
    os.makedirs(folderName)

if not os.path.exists(scratchFolderName):
    os.makedirs(scratchFolderName)

n = 18
k0 = 1
while k0 <= 16:
    k1 = 4
    while k1 <= 128:
        x = 0.9
        y = 0.9
        z = 0.9
        # print("/experiment_%d_%d_%d_%f_%f_%f.run" % (n, k0, k1, x, y, z))
        file = open(folderName + ("/experiment_%d_%d_%d_%f_%f_%f.run" % (n, k0, k1, x, y, z)), "w")
        file.write(parameters % (scratchFolderName, n, k0, k1, x, y, z))
        file.close()
        k1 *= 2
    k0 += 1
