#include <iostream>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <assert.h>

#include "bsft.h"
#include "computefourier.h"

using namespace std;

int main(int argc, char *argv[]) {
    int n = -1;
    int k0 = -1;
    int k1 = -1;
    int pow = -1;
    int numberOfTests = -1;
    double Bcst_loc = -1;
    double Bcst_est = -1;
    int loc_loops = -1;
    int est_loops = -1;
    size_t PATH_LENGTH = 128;
    char path[PATH_LENGTH + 1];
    char ch;
    while ((ch = getopt(argc, argv, "n:0:1:c:p:l:L:e:E:")) != EOF) {
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
                snprintf(path, PATH_LENGTH, "%s", optarg);
                break;
            case 'l':
                Bcst_loc = strtod(optarg, NULL);
                break;
            case 'L':
                loc_loops = strtol(optarg, NULL, 10);
                break;
            case 'e':
                Bcst_est = strtod(optarg, NULL);
                break;
            case 'E':
                est_loops = strtol(optarg, NULL, 10);
                break;
            default:
                break;
        }
    }
    if (n == -1 || k0 == -1 || k1 == -1 || pow == -1 || numberOfTests == -1) {
        cout << "Not all of the parameters are set" << endl;
        return 0;
    }
    char filename[256];
    cout << "n, k0, k1, Bcst_loc, Bcst_est, loc_loops, threshold_loops, est_loops, samples, success, time\n";
    double THRESHOLD = 0.01 * k0 * k1;

    complex_t *X = new complex_t[n];
    complex_t *X_ = new complex_t[n];
    complex_t *tmp = new complex_t[n];
    vector<double> bcsts;
//    bcsts.push_back(0.2);
//    bcsts.push_back(0.5);
    bcsts.push_back(1);
//    bcsts.push_back(2);
//    bcsts.push_back(4);
    for (loc_loops = 4; loc_loops <= 14; loc_loops += 2) {
        int threshold_loops = loc_loops - 1;
        for (est_loops = 8; est_loops <= 20; est_loops += 2) {
            for (int i_loc = 0; i_loc < bcsts.size(); i_loc++) {
                for (int i_est = 0; i_est < bcsts.size(); i_est++) {
                    Bcst_loc = bcsts[i_loc];
                    Bcst_est = bcsts[i_est];
                    double samples = 0;
                    double ticks = 0;
                    double total_succ = 0.0;
                    for (int idx = 0; idx < numberOfTests; idx++) {
                        set<int> blocks;
                        if (!load_signal_binary(path, X, blocks, n, k0, k1, idx)) {
                            continue;
                        }
                        initialize_globals(n, k1);

                        int Comb_loops = 0, W_Comb = 0;
                        int B_est, B_loc, B_thres;
                        set_sparse_fft_parameters(n,
                                                  k0 * k1,
                                                  B_est,
                                                  B_loc,
                                                  B_thres,
                                                  W_Comb,
                                                  loc_loops,
                                                  est_loops,
                                                  threshold_loops,
                                                  Comb_loops,
                                                  false,
                                                  Bcst_loc,
                                                  Bcst_est,
                                                  false);

                        clock_t prev = clock();
                        map<int, complex_t> ans = sparse_fft_wrapper(X,
                                                                     n,
                                                                     k0 * k1,
                                                                     B_est,
                                                                     B_thres,
                                                                     B_loc,
                                                                     W_Comb,
                                                                     Comb_loops,
                                                                     threshold_loops,
                                                                     loc_loops,
                                                                     est_loops);
                        ticks += clock() - prev;
                        samples += count_samples_taken(n);
                        destruct_globals(n, k1);

                        set<int> ret;
                        for (map<int, complex_t>::iterator it = ans.begin(); it != ans.end(); it++) {
                            if (cabs(it->second) < 0.001) {
                                continue;
                            }
                            ret.insert(get_block_idx(it->first, k1, n) * k1);
                        }
                        if (blocks == ret) {
                            total_succ += 1;
                        }
                    }
                    double avg_samples = samples / numberOfTests / n;
                    double avg_time = ticks / numberOfTests / CLOCKS_PER_SEC;
                    double avg_precision = total_succ / numberOfTests;
                    cout << n << ", " << k0 << ", " << k1 << ", " << Bcst_loc << ", " << Bcst_est << ", " <<
                         loc_loops << ", " << threshold_loops << ", " << est_loops << ", " << avg_samples << ", " <<
                         avg_precision << ", " << avg_time << endl;
                }
            }
        }
    }
    return 0;
}
