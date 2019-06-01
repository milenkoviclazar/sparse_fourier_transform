#!/usr/bin/python
import os

parameters = \
    """#!/bin/bash
#SBATCH --workdir /scratch/milenkov/run_heatmap_cluster_sfft/
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 40960
#SBATCH --time 24:00:00 

/home/milenkov/bsfft/heatmap_cluster_sfft -N 18 -1 4 -k %d -K %d -d 1 -l %d -L %d -p %d -P %d -b %d -B %d -v %d -V %d -n 20
    """
if not os.path.exists("sbatch_files_hikp"):
    os.makedirs("sbatch_files_hikp")


k0_lower = 2
k0_upper = 32
Bcst_loc_lower = 1
Bcst_loc_upper = 4
Bcst_est_lower = 1
Bcst_est_upper = 4
loc_loops_lower = 3
loc_loops_upper = 7
est_loops_lower = 8
est_loops_upper = 16


cnt = 0
k0 = k0_lower
while k0 <= k0_upper:
    Bcst_loc = Bcst_loc_lower
    while Bcst_loc <= Bcst_loc_upper:
        Bcst_est = Bcst_est_lower
        while Bcst_est <= Bcst_est_upper:
            for loc_loops in range(3, 8):
                for est_loops in range(8, 18, 2):
                    cnt += 1
                    file = open("sbatch_files_hikp/run_heatmap_%d_%d_%d_%d_%d.run" %
                                (k0, Bcst_loc, Bcst_est, loc_loops, est_loops), "w")
                    file.write(parameters % (k0, k0 + 1, Bcst_loc, Bcst_loc + 1, Bcst_est, Bcst_est + 1,
                                             loc_loops, loc_loops + 1, est_loops, est_loops + 1))
                    file.close()
            Bcst_est *= 2
        Bcst_loc *= 2
    k0 *= 2

print(cnt)
