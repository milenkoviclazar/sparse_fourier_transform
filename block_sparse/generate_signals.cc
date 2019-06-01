#include <iostream>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <set>
#include <cstring>
#include <cassert>

#include "bsft.h"
#include "fft.h"

using namespace std;

set<int> generate_block(complex_t X[], complex_t X_[], int n, int k0, int k1) {
    set<int> blocks;
    while (blocks.size() < k0) {
        int loc = rand_int(n / k1 - 1) + 1;
        loc *= k1;
        blocks.insert(loc);
        int lo, hi;
        get_interval_bounds(loc, k1, lo, hi, n);
        for (int i = lo; i <= hi; i++) {
            double phi = drand48() * 2 * M_PI;
            __real__ X_[i] = cos(phi);
            __imag__ X_[i] = sin(phi);
        }
    }
    fftw_dft(X, n, X_, 1);
    return blocks;
}

int main(int argc, char *argv[]) {
    int numberOfTests = 16;
    int n = 1 << 18;
    int k0 = 1;
    int k1 = 8;
    char ch;
    int pow;
    size_t PATH_LENGTH = 128;
    char path[PATH_LENGTH + 1];
    while ((ch = getopt(argc, argv, "n:0:1:c:p:")) != EOF) {
        switch (ch) {
            case 'n':
                pow = strtod(optarg, NULL);
                n = 1 << pow;
                break;
            case '0':
                k0 = strtod(optarg, NULL);
                break;
            case '1':
                k1 = strtod(optarg, NULL);
                break;
            case 'c':
                numberOfTests = strtod(optarg, NULL);
                break;
            case 'p':
                snprintf(path, PATH_LENGTH, "%s", optarg);
                break;
            default:
                break;
        }
    }
    complex_t *X = new complex_t[n];
    complex_t *X_ = new complex_t[n];
    for (int i = 0; i < numberOfTests; i++) {
        memset(X, 0, n * sizeof(complex_t));
        memset(X_, 0, n * sizeof(complex_t));
        set<int> blocks = generate_block(X, X_, n, k0, k1);
        save_signal_binary(path, X, blocks, n, k0, k1, i);
    }
    delete[] X;
    delete[] X_;
    return 0;
}
