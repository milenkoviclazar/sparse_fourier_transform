#!/usr/bin/python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sys

if(len(sys.argv) < 2):
	print("Please specify filename")
	exit(1)

filename = sys.argv[1]
print("Opening " + filename)

sns.set()

with open(filename) as f:
	lines = f.readlines()
	param = np.fromstring(lines[0], dtype=np.float64)
	data = pd.read_csv(filename, skiprows=1, header=None)
	data = data.applymap(lambda x: float(str(x).split('/')[2]))
	sns.heatmap(data)

	plt.xlabel("k1, k0 = 4, N = 2^{22}")
	plt.ylabel("p")
	
	xticks = list(map(lambda y: 2 ** y, np.linspace(3, 6+1, 4)))
	print(xticks)
	plt.xticks([0, 1, 2, 3], [8, 16, 32, 64])
	plt.yticks([0, 1, 2], [0.2, 0.15, 0.1])
plt.show()
