#!/usr/bin/python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sys

sns.set()

with open(sys.argv[1]) as f:
    lines = f.readlines()
    parameters = np.fromstring(lines[0], dtype=float, sep=',')
    n = int(parameters[0])

    sig = np.fromstring(lines[1], dtype=float, sep=',')

with open(sys.argv[1] + ".budgets") as f:
    lines = f.readlines()
    parameters = np.fromstring(lines[0], dtype=float, sep=',')
    n = int(parameters[0])
    budgets = np.fromstring(lines[1], dtype=float, sep=',')

plt.plot(sig)
plt.plot(budgets)

plt.show()
