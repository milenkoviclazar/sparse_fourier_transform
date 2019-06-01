#!/usr/bin/python
import os

parameters = \
    """#!/bin/bash
#SBATCH --workdir /scratch/milenkov/run_heatmap_cluster_bsfft_1/
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 40960
#SBATCH --time 24:00:00 

/home/milenkov/bsfft/heatmap_cluster -N %d -k %d -K %d -1 %d -n 10
    """

if not os.path.exists("sbatch_files"):
    os.makedirs("sbatch_files")

n = 18
k1 = 32
k0 = 2
while k0 <= 32:
    file = open("sbatch_files/run_heatmap_%d_%d_%d.run" % (n, k0, k1), "w")
    file.write(parameters % (n, k0, k0 + 1, k1))
    file.close()
    k0 *= 2


