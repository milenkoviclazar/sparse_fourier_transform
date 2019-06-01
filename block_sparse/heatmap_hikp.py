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

K1 = int(sys.argv[2])

max_samples = 0
best = {}
cnt = {}
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
            if k1 != K1:
                continue
            samples = float(row[8]) * n
            max_samples = max(max_samples, samples)
            precision = float(row[9])
            if (k0, samples) not in best:
                best[(samples, k0)] = precision
            else:
                best[(samples, k0)] = max(best[(samples, k0)], precision)


lim = 100
mat = np.empty((lim, 32))
mat.fill(np.NaN)

mark = [False] * len(best)
l = max_samples / lim
for i in range(lim):
    for idx, t in enumerate(best):
        k0 = int(t[1]) - 1
        if i * l <= t[0] <= (i + 1) * l:
            mark[idx] = True
            if np.isnan(mat[i][k0]):
                mat[i][k0] = best[t]
            else:
                mat[i][k0] = max(mat[i][k0], best[t])
print(mat)
print(len([idx for idx in range(len(mark)) if not mark[idx]]))
for i in range(32):
    best = 0.0
    for j in range(lim):
        if not np.isnan(mat[j][i]):
            best = max(best, mat[j][i])
        mat[j][i] = best

ax = sns.heatmap(mat)

plt.title("Best precision per samples taken")
plt.xlabel("k0")
plt.ylabel("Number of samples")
plt.xticks([i + 4 for i in range(0, 32, 4)], [i + 4 for i in range(0, 32, 4)])
plt.yticks(range(0, lim, 10), ['%.0f' % (x * l) for x in range(0, lim, 10)])
# plt.gca().yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
ax.invert_yaxis()
plt.show()