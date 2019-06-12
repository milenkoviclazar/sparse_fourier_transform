#!/usr/bin/python
import os
from datetime import datetime

timestamp_str = str(datetime.now()).replace(" ", "_").replace(":", "_")
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


folderName = "sbatch_experiment_files"
if not os.path.exists(folderName):
    os.makedirs(folderName)

if not os.path.exists(scratchFolderName):
    os.makedirs(scratchFolderName)

rng = {}
rng[4] = range(12, 33)
rng[8] = range(24, 33)
rng[16] = range(24, 33)
rng[32] = range(15, 33)
rng[64] = range(10, 33)

n = 22
k1 = 4
while k1 <= 64:
    for k0 in rng[k1]:
        x = 0.9
        y = 0.9
        z = 0.9
        # print("/experiment_%d_%d_%d_%f_%f_%f.run" % (n, k0, k1, x, y, z))
        file = open(folderName + ("/experiment_%d_%d_%d_%f_%f_%f.run" % (n, k0, k1, x, y, z)), "w")
        file.write(parameters % (scratchFolderName, n, k0, k1, x, y, z))
        file.close()
    k1 *= 2
