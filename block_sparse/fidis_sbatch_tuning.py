#!/usr/bin/python
import os

parameters = \
    """#!/bin/bash
#SBATCH --workdir /scratch/milenkov/run_tuning_18_08/
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16000
#SBATCH --time 24:00:00

/home/milenkov/fidis/bsfft/optimized/k_sparse_tuning -n %d -0 %d -1 %d -c 20 -p /scratch/milenkov/random_signals/ 
"""

folderName = "tuning_sbatch_files"
if not os.path.exists(folderName):
    os.makedirs(folderName)

if not os.path.exists("/scratch/milenkov/run_tuning_18_08/"):
    os.makedirs("/scratch/milenkov/run_tuning_18_08/")

n = 22
for k0 in range(1, 33):
    k1 = 4
    while k1 <= 64:
        file = open(folderName + ("/tune_k_sparse_%d_%d_%d.run" % (n, k0, k1)), "w")
        file.write(parameters % (n, k0, k1))
        file.close()
        k1 *= 2
