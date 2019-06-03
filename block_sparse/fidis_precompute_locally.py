#!/usr/bin/python
import os
import sys

if len(sys.argv) < 3:
    print("expected arguments: signalExecPath, signalOutPath")
    quit()

signalExecPath = sys.argv[1] + "/"
signalOutPath = sys.argv[1] + "/"
if not os.path.exists(signalOutPath):
    os.makedirs(signalOutPath)
for n in range(22, 23):
    for k0 in range(1, 33):
        k1 = 2
        while k1 <= 64:
            print(n, k0, k1)
            if k0 * k1 < pow(2, n) / 2:
                os.system("%s -n %d -0 %d -1 %d -c 20 -p %s" % (signalExecPath, n, k0, k1, signalOutPath))
            k1 *= 2
