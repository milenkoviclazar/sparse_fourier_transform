import csv
import sys
from os import listdir
from os.path import isfile, join
import math
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

path = sys.argv[1] + '/'
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

N = int(pow(2, int(sys.argv[2])))
K0 = int(sys.argv[3])
# 0, 1 , 2 ,     3   ,    4    ,     5    ,        6       ,      7   ,   8    ,   9    ,  10
# n, k0, k1, Bcst_loc, Bcst_est, loc_loops, threshold_loops, est_loops, samples, success, time
max_samples = 0
best = {}
best_parameters = {}
for file in onlyfiles:
    if file[0] == ".":
        continue
    with open(path + file, 'r') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            if len(row) != 11:
                # print(row)
                continue
            try:
                n = int(row[0])
                k0 = int(row[1])
                k1 = int(row[2])
            except ValueError:
                # print(row)
                continue
            if n != N:
                continue
            if k0 == K0 and float(row[9]) > 0.9:
                print(row)

            samples = float(row[8]) * n
            max_samples = max(max_samples, samples)
            precision = float(row[9])
            if ((k0, samples) not in best) or (precision > best[(samples, k0)]):
                best[(samples, k0)] = precision
                best_parameters[(samples, k0)] = [float(row[i]) for i in range(3, 8)]

exit()

COLS = 32
lim = 100
mat = np.empty((lim, COLS))
mat.fill(np.NaN)

param_idx = 0
mark = [False] * len(best)
l = max_samples / lim
for i in range(lim):
    for idx, t in enumerate(best):
        k0 = int(t[1]) - 1
        if i * l <= t[0] <= (i + 1) * l:
            mark[idx] = True
            if np.isnan(mat[i][k0]):
                # mat[i][k0] = best_parameters[t][param_idx] # best parameters
                mat[i][k0] = best[t]
            else:
                # mat[i][k0] = min(mat[i][k0], best_parameters[t][param_idx])
                mat[i][k0] = max(mat[i][k0], best[t])

print(mat)
print(len([idx for idx in range(len(mark)) if not mark[idx]]))
for i in range(COLS):
    best = 0.0
    for j in range(lim):
        if not np.isnan(mat[j][i]):
            best = max(best, mat[j][i])
        # mat[j][i] = best

ax = sns.heatmap(mat)

plt.title("Best precision per samples taken")
plt.xlabel("k0")
plt.ylabel("Number of samples")
SPACING = 4
plt.xticks([i + SPACING for i in range(0, COLS, SPACING)], [i + SPACING for i in range(0, COLS, SPACING)])
plt.yticks(range(0, lim, 10), ['%.0f' % (x * l) for x in range(0, lim, 10)])
# plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
ax.invert_yaxis()
plt.show()