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
	parameters = np.fromstring(lines[0], dtype=float, sep=',')
	n = int(parameters[0])

	x_real = np.fromstring(lines[1], dtype=float, sep=',')
	x_imag = np.fromstring(lines[2], dtype=float, sep=',')
	x_abs2 = np.array(list(map(lambda x: x[0] ** 2 + x[1] ** 2, list(zip(x_real, x_imag)))))

	fig, axes = plt.subplots(3, 1)
	axes[0].plot(x_real)
	axes[0].set_title("real(f)")
	axes[1].plot(x_imag)
	axes[1].set_title("imag(f)")
	axes[2].plot(x_abs2)
	axes[2].set_title("abs2(f)")

	plt.show()
