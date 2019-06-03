#!/usr/bin/python
import os
import sys

if len(sys.argv) < 5:
    print("expected arguments: filterExecPath, filterOutPath, signalExecPath, signalOutPath")
    quit()

filterExecPath = sys.argv[1]
filterOutPath = sys.argv[2] + "/"
if not os.path.exists(filterOutPath):
    os.makedirs(filterOutPath)

for n in range(18, 19):
    b = 2
    while b <= pow(2, n) / 2:
        os.system("%s -n %d -b %d -p %s" % (filterExecPath, n, b, filterOutPath))
        b *= 2

signalExecPath = sys.argv[3]
signalOutPath = sys.argv[4] + "/"

if not os.path.exists(signalOutPath):
    os.makedirs(signalOutPath)

for n in range(18, 19):
    k0 = 1
    while k0 <= 32:
        k1 = 2
        while k1 <= 64:
            print(n, k0, k1)
            if k0 * k1 < pow(2, n) / 2:
                os.system("%s -n %d -0 %d -1 %d -c 20 -p %s" % (signalExecPath, n, k0, k1, signalOutPath))
            k1 *= 2
        k0 += 1
