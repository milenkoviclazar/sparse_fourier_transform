#include <iostream>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "bsft.h"
#include "computefourier.h"

using namespace std;

int main(int argc, char *argv[]) {
    int n = -1;
    int k0 = -1;
    int k1 = -1;
    int pow = -1;
    int numberOfTests = 1;

    // TODO: Should be removed. It is a kind of a heuristic and should not work for general signals.
    double hash_percentage = 0.9;
    double sample_percentage = 0.9;
    double val_percentage = 0.9;

    size_t PATH_LENGTH = 128;
    char filterPath[PATH_LENGTH + 1];
    memset(filterPath, 0, sizeof(filterPath));
    char signalPath[PATH_LENGTH + 1];
    memset(signalPath, 0, sizeof(signalPath));
    char ch;

    while ((ch = getopt(argc, argv, "n:0:1:c:p:q:x:y:z:")) != EOF) {
        switch (ch) {
            case 'n':
                pow = strtol(optarg, NULL, 10);
                n = 1 << pow;
                break;
            case '0':
                k0 = strtol(optarg, NULL, 10);
                break;
            case '1':
                k1 = strtol(optarg, NULL, 10);
                break;
            case 'c':
                numberOfTests = strtol(optarg, NULL, 10);
                break;
            case 'p':
                snprintf(signalPath, PATH_LENGTH, "%s", optarg);
                break;
            case 'q':
                snprintf(filterPath, PATH_LENGTH, "%s", optarg);
                break;
            case 'x':
                hash_percentage = strtof(optarg, NULL);
                break;
            case 'y':
                sample_percentage = strtof(optarg, NULL);
                break;
            case 'z':
                val_percentage = strtof(optarg, NULL);
                break;
            default:
                break;
        }
    }
    if (n == -1 || k0 == -1 || k1 == -1 || pow == -1 || filterPath[0] == 0 || signalPath[0] == 0
        || numberOfTests == -1) {
        cout << "Not all of the parameters are set" << endl;
        return 0;
    }
    complex_t *X = new complex_t[n];
    complex_t *X_ = new complex_t[n];
    double theta_val = 0.01;
    int F_val = 10;
    int F_sample = 10;
    int F_loc = 10;
#ifdef __APPLE__
    int min_samples = n;
#endif
    cout << "n, k0, k1, B_loc, B_val, iter_loc, iter_budget, iter_val, hash_p, val_p, sample_p, "
            "avg_samples, avg_time, succ_prob" << endl;

    for (int B_loc_mult = 1; B_loc_mult < 8; B_loc_mult++) {
        int B_loc = round_to_power2(B_loc_mult * k0);
        for (int B_val_mult = 1; B_val_mult < 128; B_val_mult += 2) {
            int B_val = round_to_power2(B_val_mult * k0 * k1);

#ifdef __APPLE__
            Filter H_hash(n / k1, B_loc, F_loc);
            Filter H_sample(n, n / k1, F_sample);
            Filter G_val(n, B_val, F_val);
#else
            Filter H_hash(filterPath, n / k1, B_loc, F_loc);
            H_hash.size = min(n / 2, 8 * B_loc * F_loc);
            Filter H_sample(filterPath, n, n / k1, F_sample);
            H_sample.size = min(n / 2, 8 * n / k1 * F_sample);
            Filter G_val(filterPath, n, B_val, F_val);
            G_val.size = min(n / 2, 8 * B_val * F_val);
            if (H_hash.n == -1 || H_sample.n == -1 || G_val.n == -1) {
                cerr << "Filter problem" << endl;
                return 0;
            }
#endif

//            for (double hash_percentage = 0.8; hash_percentage <= 1.0; hash_percentage += 0.5) {
//
                H_hash.resize(hash_percentage);

//
//                for (double sample_percentage = 0.6; sample_percentage <= 1.0; sample_percentage += 0.5) {
                    H_sample.resize(sample_percentage);
//
//                    for (double val_percentage = 0.75; val_percentage <= 1.0; val_percentage += 0.5) {
                        G_val.resize(val_percentage);

                        int iter_budget = k0;
                        for (int iter_loc = 0; iter_loc < 10; iter_loc += 3) {
                            for (int iter_val = 1; iter_val < 40; iter_val += 3) {
                                double ticks = 0;
                                double succ = 0;
                                double samples = 0;
                                for (int i = 0; i < numberOfTests; i++) {
                                    set<int> trueBlocks;
                                    if (!load_signal_binary(signalPath, X, trueBlocks, n, k0, k1, i)) {
                                        continue;
                                    }

                                    initialize_globals(n, k1);
                                    set_parameters(n, k0, k1);
                                    clock_t prev = clock();
                                    set<int> guessBlocks;
                                    guessBlocks = bsft_location(X, n, k0, k1, G_val, H_sample, H_hash, B_loc, B_val,
                                                                iter_loc, iter_budget, iter_val, theta_val);
                                    ticks += clock() - prev;
                                    samples += count_samples_taken(n);
                                    succ += (guessBlocks == trueBlocks);
                                    destruct_globals(n, k1);
                                }
                                double avg_time = ticks / numberOfTests / CLOCKS_PER_SEC;
                                double succ_prob = succ / numberOfTests;
                                double avg_samples = samples / numberOfTests;
#ifdef __APPLE__
//                    if (succ_prob == 1.0) {
//                        min_samples = min(min_samples, (int) avg_samples);
//                    }
//                    cout << min_samples << endl;
#endif
                                cout << n << ", " << k0 << ", " << k1 << ", " << B_loc << ", " << B_val
                                     << ", " <<
                                     iter_loc << ", " << iter_budget << ", " << iter_val << ", " <<
                                     hash_percentage << ", " << val_percentage
                                     <<
                                     ", " << sample_percentage << ", " << avg_samples << ", " << avg_time << ", "
                                     << succ_prob << endl;
                            }
                        }
//                    }
//                }
//            }

            H_hash.clr();
            H_sample.clr();
            G_val.clr();
        }
    }
    delete[] X;
    delete[] X_;
    return 0;
}
