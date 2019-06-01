#include <iostream>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "bsft.h"

using namespace std;

int main(int argc, char *argv[]) {
    int pow = -1;
    int n = -1;
    int B = -1;
    int F = 10;
    char ch;
    size_t PATH_LENGTH = 128;
    char filterPath[PATH_LENGTH + 1];
    memset(filterPath, 0, sizeof(filterPath));
    while ((ch = getopt(argc, argv, "n:b:f:p:")) != EOF) {
        switch (ch) {
            case 'n':
                pow = strtod(optarg, NULL);
                n = 1 << pow;
                break;
            case 'b':
                B = strtod(optarg, NULL);
                break;
            case 'f':
                F = strtol(optarg, NULL, 10);
            case 'p':
                snprintf(filterPath, PATH_LENGTH, "%s", optarg);
                break;
            default:
                break;
        }
    }
    if (n == -1 || B == -1 || filterPath[0] == 0) {
        cerr << "Not all the arguments are set";
        return 0;
    }
    Filter G(n, B, F);
    G.dump_binary(filterPath);
    return 0;
}
