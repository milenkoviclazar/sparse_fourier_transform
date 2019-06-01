#!/usr/bin/python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sys

if(len(sys.argv) < 2):
	print("Please specify filename")
	exit(1)

nb = 0
filename = sys.argv[1]
print("Opening " + filename)

if(len(sys.argv) > 2):
	nb = int(sys.argv[2])

sns.set()

with open(filename) as f:
	lines = f.readlines()
	parameters = np.fromstring(lines[0], dtype=float, sep=',')
	n = parameters[0]
	k0 = parameters[1]
	k1 = parameters[2]

	x_real = np.fromstring(lines[1], dtype=float, sep=',')
	x_imag = np.fromstring(lines[2], dtype=float, sep=',')
	x_abs2 = np.fromstring(lines[3], dtype=float, sep=',')
	energies = np.fromstring(lines[4], dtype=float, sep=',')
	budgets = np.fromstring(lines[5], dtype=float, sep=',')

	fig, axes = plt.subplots(3, 2)
	axes[0, 0].plot(x_real)
	axes[0, 0].set_title("real(X)")
	axes[1, 0].plot(x_imag)
	axes[1, 0].set_title("imag(X)")
	axes[2, 0].plot(x_abs2)
	axes[2, 0].set_title("abs2(X)")
	axes[0, 1].stem(energies)
	axes[0, 1].set_title("energy estimates")
	axes[1, 1].stem(budgets)
	axes[1, 1].set_title("Budgets allocations")

	axes[2, 1].plot(x_abs2)
	if nb != 0:
		sorted_budgets = np.argsort(-budgets)
		m = 0
		_, counts = np.unique(budgets, return_counts=True)
		np.sort(counts)
		for j in range(nb):
			m += counts[-(j+1)]	
		maxi = np.max(x_abs2)
		s=0
		total=np.sum(x_abs2)
		for i in sorted_budgets[:m]:
			i1 = int(i*n/(2*k1))
			i2 = int((i+1)*n/(2*k1))
			for j in range(i1, i2):
				s+=x_abs2[j]
			axes[2, 1].fill([i1, i2, i2, i1], [0, 0, maxi, maxi], alpha=.5, color='r')

		axes[2, 1].text(0, 500, np.around(s/total, 2))
	
	plt.show()
