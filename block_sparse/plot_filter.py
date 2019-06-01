#!/usr/bin/python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sys

if len(sys.argv) < 2:
    print("Please specify filename")
    exit(1)

filename = sys.argv[1]
print("Opening " + filename)

sns.set()

with open(filename) as f:
    lines = f.readlines()
    parameters = np.fromstring(lines[0], dtype=int, sep=',')
    # print(parameters)
    n = parameters[0]
    # B = parameters[1]
    # size = parameters[2]

    x_real = np.fromstring(lines[1], dtype=float, sep=',')
    x_imag = np.fromstring(lines[2], dtype=float, sep=',')
    x_abs2 = np.copy(x_real)
    for i in range(len(x_real)):
        x_abs2[i] = x_real[i] * x_real[i] + x_imag[i] * x_imag[i]

    X_real = np.fromstring(lines[3], dtype=float, sep=',')
    X_imag = np.fromstring(lines[4], dtype=float, sep=',')
    X_abs2 = np.copy(X_imag)
    for i in range(len(x_real)):
        X_abs2[i] = X_real[i] * X_real[i] + X_imag[i] * X_imag[i]

    fig, axes = plt.subplots(3, 2)
    axes[0, 0].plot(x_real)

    axes[0, 0].set_title("Re(time)")
    axes[1, 0].plot(x_imag)
    axes[1, 0].set_title("Im(time)")
    axes[2, 0].plot(x_abs2)
    axes[2, 0].set_title("abs2(time)")
    axes[0, 1].plot(X_real)
    axes[0, 1].set_title("Re(freq)")

    axes[1, 1].plot(X_imag)
    axes[1, 1].set_title("Im(freq)")
    axes[2, 1].plot(X_abs2)
    axes[2, 1].set_title("abs2(freq")

    plt.show()
