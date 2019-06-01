#!/usr/bin/python
import os

parameters = \
    """#!/bin/bash
#SBATCH --workdir /scratch/milenkov/run_experiments_23_01_2019/
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 50GB
#SBATCH --time 06:00:00

/home/milenkov/deneb/bsfft/optimized/run_experiment -n %d -0 %d -1 %d -c 20 -x %f -y %f -z %f \
-p /scratch/milenkov/random_signals/ \
-q /scratch/milenkov/filters10/
"""


folderName = "experiment_sbatch_files"
if not os.path.exists(folderName):
    os.makedirs(folderName)


n = 18
k0 = 2
while k0 <= 32:
    k1 = 4
    while k1 <= 64:
        x = 0.6
        while x <= 1.0:
            y = 0.6
            while y <= 1.0:
                z = 0.6
                while z <= 1.0:
                    # print("/experiment_%d_%d_%d_%f_%f_%f.run" % (n, k0, k1, x, y, z))
                    file = open(folderName + ("/experiment_%d_%d_%d_%f_%f_%f.run" % (n, k0, k1, x, y, z)), "w")
                    file.write(parameters % (n, k0, k1, x, y, z))
                    file.close()
                    z += 0.05
                y += 0.05
            x += 0.05
        k1 *= 2
    k0 += 1
