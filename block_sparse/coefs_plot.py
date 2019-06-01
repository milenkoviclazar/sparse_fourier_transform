#!/usr/bin/python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sys

if(len(sys.argv) < 3):
	print("Please specify filename")
	exit(1)

filename = sys.argv[1]
filename2 = sys.argv[2]
print("Opening " + filename)

sns.set()
fig, axes = plt.subplots(5, 1)

with open(filename) as f:
	lines = f.readlines()
	parameters = np.fromstring(lines[0], dtype=float, sep=',')
	n = parameters[0]

	x_real = np.fromstring(lines[1], dtype=float, sep=',')
	x_imag = np.fromstring(lines[2], dtype=float, sep=',')
	x_abs2 = np.array(list(map(lambda x: x[0] ** 2 + x[1] ** 2, list(zip(x_real, x_imag)))))
#	x_abs2 = np.fromstring(lines[3], dtype=float, sep=',')

	axes[0].plot(x_real)
	axes[0].set_title("real(f)")
	axes[1].plot(x_imag)
	axes[1].set_title("imag(f)")
	axes[2].plot(x_abs2)
	axes[2].set_title("abs2(f)")

print("Opening " + filename2)
df = pd.read_csv(filename2, header=None)
x = np.transpose(df.values)
print(x)
print(np.min(np.diff(x[0])))
axes[3].plot(x_abs2)
axes[3].stem(x[0], x[1])
if(len(sys.argv) > 3):
	df = pd.read_csv(sys.argv[3])
	x = np.transpose(df.values)
	x_abs2 = np.array(list(map(lambda x: x[0] ** 2 + x[1] ** 2, list(zip(x_real, x_imag)))))
	axes[4].plot(x_abs2)
	axes[4].plot(x[0], x[1])

plt.show()
