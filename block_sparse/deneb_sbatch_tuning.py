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

/home/milenkov/deneb/sparse_fourier_transform/block_sparse/k_sparse_tuning -n %d -0 %d -1 %d -c 20 \
-p /scratch/milenkov/signals/"""


folderName = "sbatch_tuning_files"
if not os.path.exists(folderName):
    os.makedirs(folderName)

if not os.path.exists(scratchFolderName):
    os.makedirs(scratchFolderName)

n = 22
for k0 in range(1, 33):
    k1 = 4
    while k1 <= 64:
        file = open(folderName + ("/tune_k_sparse_%d_%d_%d.run" % (n, k0, k1)), "w")
        file.write(parameters % (scratchFolderName, n, k0, k1))
        file.close()
        k1 *= 2