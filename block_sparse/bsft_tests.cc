#include <set>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <fstream>

#include "bsft.h"
#include "computefourier.h"

typedef double _Complex complex_t;

using namespace std;

string BSFT_TERMINAL = "\033[33m";
string HIKP_TERMINAL = "\033[36m";
string ESC_TERMINAL = "\033[0m";

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

set<int> generate_location_testcase(complex_t X[], complex_t X_[], int n, int k0, int k1, int peaks) {
    set<int> blocks;
    set<int> energies;
    while (blocks.size() < k0) {
        int loc = rand_int(n / k1 - 1) + 1;
        loc *= k1;
        blocks.insert(loc);
    }
    while (energies.size() < peaks) {
        energies.insert(rand_int(n));
    }
    complex_t omega = cexp(2 * I * M_PI / n);
    for (set<int>::iterator bl = blocks.begin(); bl != blocks.end(); bl++) {
        int lo, hi;
        get_interval_bounds(*bl, k1, lo, hi, n);
        for (int i = lo; i <= hi; i++) {
            for (set<int>::iterator en = energies.begin(); en != energies.end(); en++) {
                X_[i] += cpow(omega, n - *en * (i % n));
            }

        }
    }
    fftw_dft(X, n, X_, 1);
    return blocks;
}

set<int> indices_slow;

void slow_hash(
        complex_t *out,
        complex_t *X,
        int n,
        long long B,
        Filter G,
        long long sigma,
        long long delta) {
    for (long long j = 0; j < B; j++) {
        out[j] = 0;
        for (size_t idx = 0; idx < n / B; idx++) {
            indices_slow.insert((sigma * (delta + j + B * idx)) % n);
            out[j] += X[(sigma * (delta + j + B * idx)) % n] * G.time[(j + B * idx) % n];
        }
        out[j] *= (double) B / (double) n;
    }
}

void slow_hash_to_bins(
        int n,
        complex_t *X,
        complex_t *U,
        complex_t *U_,
        Filter &G,
        int B,
        int sigma,
        int Delta) {
    slow_hash(U, X, n, B, G, sigma, Delta);
    fftw_dft(U_, B, U);
    for (int i = 0; i < B; i++) {
        U_[i] /= B;
    }
}

double slow_hash_to_bins_reduced(
        int n,
        complex_t *Z,
        complex_t *U,
        complex_t *U_,
        Filter &G,
        int B,
        int sigma,
        int delta) {
    slow_hash(U, Z, n, B, G, sigma, delta);
    fftw_dft(U_, B, U);
    double ret = 0.0;
    for (int i = 0; i < B; i++) {
        ret += cabs2(U_[i] / B);
    }
    return ret;
}

double test_h2b(int n, int k0, int k1) {
    initialize_globals(n, k1);
    complex_t *X = new complex_t[n];
    memset(X, 0, n * sizeof(complex_t));
    complex_t *X_ = new complex_t[n];
    memset(X_, 0, n * sizeof(complex_t));
    set<int> blocks = generate_block(X, X_, n, k0, k1);

    double max_error = 0.0;
    for (int B = 2; B < n; B <<= 1) {
        cout << "B: " << B << endl;
        Filter G(n, B, 10);
        cout << "Relative size: " << 1.0 * G.size / n << endl;
        double filter_energy = 0.0;
        double filter_ignored = 0.0;
        for (int i = 0; i < n; i++) {
            if (i < G.size || i + G.size > n) {
                filter_energy += cabs2(G.time[i]);
            } else {
                filter_ignored += cabs2(G.time[i]);
            }
        }
        cout << "ignored part: " << filter_ignored / (filter_ignored + filter_energy) << endl;
        complex_t *normalU = new complex_t[B];
        complex_t *normalU_ = new complex_t[B];
        complex_t *slowU = new complex_t[B];
        complex_t *slowU_ = new complex_t[B];
        for (int iter = 0; iter < 10; iter++) {
            int sigma = rand_int(n / k1 / 2) * 2 + 1;
            int Delta = rand_int(n / k1);
            hash_to_bins(n, X, normalU, normalU_, G, B, sigma, Delta);
            slow_hash_to_bins(n, X, slowU, slowU_, G, B, sigma, Delta);
            double error = 0.0;
            for (int i = 0; i < B; i++) {
                error += cabs2(normalU_[i] - slowU_[i]);
            }
            error /= B;
            cout << "The average error per bin is " << error / B << endl;
            max_error = max(max_error, error);
        }
    }
    delete[] X;
    delete[] X_;
    destruct_globals(n, k1);
    return max_error;
}

double test_h2b_shape(int n, int k0, int k1) {
    initialize_globals(n, k1);
    srand(time(0));
    complex_t *X = new complex_t[n];
    memset(X, 0, n * sizeof(complex_t));
    complex_t *X_ = new complex_t[n];
    memset(X_, 0, n * sizeof(complex_t));
    set<int> blocks = generate_block(X, X_, n, k0, k1);
    int B = round_to_power2(k0 * k1 * 16);
    complex_t *U = new complex_t[B];
    memset(U, 0, B * sizeof(complex_t));
    complex_t *U_ = new complex_t[B];
    memset(U_, 0, B * sizeof(complex_t));
    int F = 10;
    Filter G(n, B, F);
    int sigma = rand_int(n / k1 / 2) * 2 + 1;
    int Delta = rand_int(n / k1);
    slow_hash_to_bins(n, X, U, U_, G, B, sigma, Delta);
    save_signal("OUT/h2b", U_, B);
    destruct_globals(n, k1);
    delete[] X;
    delete[] X_;
    delete[] U;
    delete[] U_;
    return 0;
}

double test_reduced_h2b(int n, int k0, int k1) {
    initialize_globals(n, k1);
    cout << n << endl;
    complex_t *X = new complex_t[n];
    memset(X, 0, n * sizeof(complex_t));
    complex_t *X_ = new complex_t[n];
    memset(X_, 0, n * sizeof(complex_t));
    set<int> blocks = generate_block(X, X_, n, k0, k1);

    int F = 10;
    Filter G_downsample(n, n / k1, F);
//    G_downsample.dump("OUT/filtercic");
//    exit(0);
    complex_t **Z = new complex_t *[2 * k1];
    for (int i = 0; i < 2 * k1; i++) {
        Z[i] = new complex_t[n / k1];
        for (int j = 0; j < n / k1; j++) {
            Z[i][j] = get_downsampled(X, G_downsample, i, j, k1, n);
        }
    }

    for (int B = 2; B < n / k1; B <<= 1) {
        Filter G_hash(n / k1, B, F);
        complex_t *U = new complex_t[B];
        complex_t *U_ = new complex_t[B];
        complex_t *slowU = new complex_t[B];
        complex_t *slowU_ = new complex_t[B];
        double error = 0.0;
        for (int r = 0; r < 2 * k1; r++) {
            int sigma = 1;
            int Delta = 1;
            double fast = hash_to_bins_reduced(n, k1, r, X, U, U_, G_hash, G_downsample, B, sigma, Delta);
            double slow = slow_hash_to_bins_reduced(n / k1, Z[r], slowU, slowU_, G_hash, B, sigma, Delta);
        }
        cout << B << " " << error << endl;
        G_hash.clr();
        delete[] U;
        delete[] U_;
        delete[] slowU;
        delete[] slowU_;
    }

    G_downsample.clr();
    delete[] X;
    delete[] X_;
    for (int i = 0; i < 2 * k1; i++) {
        delete[] Z[i];
    }
    delete[] Z;
    destruct_globals(n, k1);
    return 0.0;
}

double test_location(int n, int k0, int k1, double delta, double p) {
    srand(time(0));
    initialize_globals(n, k1);
    complex_t *X = new complex_t[n];
    memset(X, 0, n * sizeof(complex_t));
    complex_t *X_ = new complex_t[n];
    memset(X_, 0, sizeof(complex_t));
//    set<int> locations = generate_block(X, X_, n, k0, k1);
    X_[0] = 1;
    fftw_dft(X, n, X_, 1);
    save_signal("OUT/X", X, n);
    save_signal("OUT/X_", X_, n);

    int B = k0 << 6;
    Filter H_hash(n, B, 10);
    Filter H_downsample(n, n / k1, 10);
    int loc_iter = 10;
    int budget_iter = 10;
    set<int> result;
//    multi_block_locate(X, H_hash, H_downsample, result, n, k0, k1, delta, p, B, loc_iter, budget_iter);

    delete[] X;
    delete[] X_;
    destruct_globals(n, k1);
    return 0;
}


complex_t get_downsampled_old(complex_t *X, Filter &G, int r, int j, int k1, int n) {
//    if (*(unsigned long long *) (&cache[r][j]) != 0xfefefefefefefefe) {
//        return cache[r][j];
//    }
    int nk1 = n / k1;
    complex_t result = 0.0;
    for (int i = 0; i < k1; i++) {
        int idxG = j + nk1 * i;
        while (idxG >= n) {
            idxG -= n;
        }
        int idxX = idxG + nk1 / 2 * r;
        while (idxX >= n) {
            idxX -= n;
        }
        result += G.time[idxG] * take_sample(X, idxX);
    }
//    cache[r][j] = result / k1;
//    return cache[r][j];
    return 0;
}

int cnt = 0;

complex_t get_downsampled_new1(complex_t *X, Filter &G, int r, int j, int k1, int n) {
    int nk1 = n / k1;
    complex_t result = 0.0;
    for (int i = 0; i < k1; i++) {
        int idxG = j + nk1 * i;
        while (idxG >= n) {
            idxG -= n;
        }
        int idxX = idxG + nk1 / 2 * r;
        while (idxX >= n) {
            idxX -= n;
        }
        if (idxG < G.size || idxG + G.size >= n) {
            cnt++;
            result += G.time[idxG] * X[idxX]; // TODO: take sample
        }
    }
    return result / k1;
}

void test_downsampling(int n, int k0, int k1, double delta, double p) {
    complex_t *X = new complex_t[n];
    memset(X, 0, n * sizeof(complex_t));
    complex_t *X_ = new complex_t[n];
    memset(X_, 0, sizeof(complex_t));
//    set<int> locations = generate_block(X, X_, n, k0, k1);
    for (int i = 0; i < n; i++) {
        X[i] = 1.0;
    }
    Filter G(n, n / k1, 10);

    cout << "SAMPLES: " << cnt / (2 * n) << endl;

    double maxi = 0.0;
    double total = 0.0;
    complex_t **d1 = new complex_t *[2 * k1];
    complex_t **d2 = new complex_t *[2 * k1];
    for (int i = 0; i < 2 * k1; i++) {
        d1[i] = new complex_t[n / k1];
        d2[i] = new complex_t[n / k1];
    }
    initialize_globals(n, k1);
    clock_t prev = clock();
    for (int r = 0; r < 2 * k1; r++) {
        for (int j = 0; j < n / k1; j++) {
            d1[r][j] = get_downsampled(X, G, r, j, k1, n);
        }
    }
    cout << fixed << setprecision(5) << "downsampled previous: " << double(clock() - prev) / CLOCKS_PER_SEC << endl;
    destruct_globals(n, k1);

//    int cnt = 0;
    prev = clock();
    for (int r = 0; r < 2 * k1; r++) {
        for (int j = 0; j < n / k1; j++) {
            d2[r][j] = get_downsampled_new1(X, G, r, j, k1, n);
        }
    }
//    cout << 1.0 * cnt / (n * 2 * k1) << endl;


    cout << fixed << setprecision(5) << "downsampled new: " << double(clock() - prev) / CLOCKS_PER_SEC << endl;
    for (int r = 0; r < 2 * k1; r++) {
        for (int j = 0; j < n / k1; j++) {
            double diff = cabs2(d1[r][j] - d2[r][j]) / cabs2(d1[r][j]);
            maxi = max(maxi, diff);
            total += diff;
        }
    }
    cout << "MAX diff: " << maxi << endl;
    cout << "AVG diff: " << total / n << endl;
}

void test_binary() {
    int n = 1 << 22;
    long long *arr = new long long[n];
    for (int i = 0; i < n; i++) {
        arr[i] = i;
    }
    size_t size = n * sizeof(long long);
    char *buff = new char[size];
    memcpy(buff, arr, size);
    ofstream ofs("OUT/binary", ofstream::binary);
    ofs.write(buff, size);
    ofs.close();
    delete[] buff;

    ifstream ifs("OUT/binary", ifstream::binary);
    buff = new char[size];
    memset(buff, 0, size);
    ifs.read(buff, size);
    memset(arr, 0, size);
    memcpy(arr, buff, size);
//    for (int i = 0; i < n; i++) {
//        cout << arr[i] << " ";
//    }
}


void minimize_samples(int n, int k0, int k1) {
    complex_t *X = new complex_t[n];
    memset(X, 0, n * sizeof(complex_t));
    complex_t *X_ = new complex_t[n];
    memset(X_, 0, n * sizeof(complex_t));
    set<int> blocks;
    int B_loc = 2 * k0;
    int B_val = 256 * k0 * k1;
    Filter H_hash(n / k1, B_loc, 10);
    Filter H_sample(n, n / k1, 10);
    Filter G_val(n, B_val, 10);
    double sampling_filter_energy = 0.0;
    for (int i = 0; i < n / 2; i++) {
        sampling_filter_energy += cabs(H_sample.time[i]);
    }
    double energy_sum = 0.0;
    for (H_sample.size = 0; H_sample.size < n / 2; H_sample.size++) {
        energy_sum += cabs(H_sample.time[H_sample.size]);
        if (energy_sum >= sampling_filter_energy * 0.9) {
            break;
        }
    }

    int iter_loc = 3;
    int iter_budget = 2;
    int iter_val = 6;
    double theta_val = 1;
    int NUMBER_OF_TESTS = 5;
    int cnt = 0;
    for (int i = 0; i < NUMBER_OF_TESTS; i++) {
        blocks = generate_block(X, X_, n, k0, k1);
        map<int, complex_t> result;
        initialize_globals(n, k1);
        set_parameters(n, k0, k1);
        result = bsft(X, n, k0, k1, G_val, H_sample, H_hash, B_loc, B_val, iter_loc, iter_budget, iter_val, theta_val);
        destruct_globals(n, k1);
        set<int> recognized;
        complex_t *tmp = new complex_t[n];
        memset(tmp, 0, n * sizeof(complex_t));
        for (map<int, complex_t>::iterator it = result.begin(); it != result.end(); it++) {
//            cout << it->first << " ";
            tmp[it->first] = it->second;
            recognized.insert(get_block_idx(it->first, k1, n) * k1);
        }
        save_signal("OUT/tmp", tmp, n);
        cout << "PLOT WAIT..";
        getchar();
        for (set<int>::iterator it = recognized.begin(); it != recognized.end(); it++) {
            cout << *it << " ";
        }
        cout << endl;
        for (set<int>::iterator it = blocks.begin(); it != blocks.end(); it++) {
            cout << *it << " ";
        }
        cout << endl << endl;
//        cout << blocks.size() << " " << recognized.size() << " ";
//        cout << (recognized == blocks ? "YES": "NO") << endl;

        cnt += (recognized == blocks);
    }
    cout << 1.0 * cnt / NUMBER_OF_TESTS << endl;

}

void tune_pruning(int n, int k0, int k1) {
    complex_t *X = new complex_t[n];
    memset(X, 0, n * sizeof(complex_t));
    complex_t *X_ = new complex_t[n];
    memset(X_, 0, n * sizeof(complex_t));
    int numberOfSimulations = 40;
    int min_samples = n;
    for (int B_val_mult = 32; B_val_mult < 128; B_val_mult <<= 1) {
        int B_val = B_val_mult * k0 * k1;
        Filter G_val(n, B_val, 2);
//        G_val.dump("OUT/filter");
//        cout << "WAIT PLOT...." << endl;
//        getchar();
        for (double filter_percent = 0.7; filter_percent < 1.0; filter_percent += 0.05) {
            double total_filter_energy = 0.0;
            for (int i = 0; i < n / 2; i++) {
                total_filter_energy += cabs(G_val.time[i]);
            }
            double curr_sum = 0.0;
            for (G_val.size = 0; G_val.size < n / 2; G_val.size++) {
                curr_sum += cabs(G_val.time[G_val.size]);
                if (curr_sum > filter_percent * total_filter_energy) {
                    break;
                }
            }
            for (int iter_val = 1; iter_val < 20; iter_val += 2) {
                double samples = 0;
                double cnt = 0;
                for (int i = 0; i < numberOfSimulations; i++) {
                    memset(X, 0, n * sizeof(complex_t));
                    memset(X_, 0, n * sizeof(complex_t));
                    set<int> blocks = generate_block(X, X_, n, k0, k1);
                    set<int> candidates(blocks.begin(), blocks.end());
                    for (int j = k1; j < n; j += k1) {
                        if (drand48() < 0.001) {
                            candidates.insert(j);
                        }
                    }
                    map<int, complex_t> result;
                    initialize_globals(n, k1);
                    set<int> ans = prune_location(X, n, k0, k1, candidates, G_val, B_val, iter_val, 0.01);
//                    cout << blocks.size() << " " << ans.size() << endl;
//                    for (set<int>::iterator it = ans.begin(); it != ans.end(); it++) {
//                        cout << *it << " "
//                    }
                    cnt += (blocks == ans);
                    samples += count_samples_taken(n);
                    destruct_globals(n, k1);
                }
                double avg_samples = samples / numberOfSimulations;
                double avg_precision = cnt / numberOfSimulations;
                if (avg_precision >= 0.95) {
                    min_samples = min(min_samples, (int)avg_samples);
                }
                cout << iter_val << ", " << B_val_mult << ", " << filter_percent << ", " <<
                     avg_precision << ", " << avg_samples << endl;
            }
        }
    }
    cout << min_samples << endl;
}

int main() {
    int n = 255;

    double Bcst_loc = 1.0;
    int k = 100;

    real_t BB_loc = (unsigned) (Bcst_loc * sqrt((double) n * k / (log2(n))));
    double lobefrac_loc = 0.5 / (BB_loc);
    int w_loc;
    int b_loc = int(1.2 * 1.1 * ((double) n / BB_loc));

    lobefrac_loc = 0.11;
    double tolerance_loc = 2e-9;

    complex_t *filtert = make_dolphchebyshev_t(lobefrac_loc, tolerance_loc, w_loc);
    FilterHikp f = make_multiple_t(filtert, w_loc, n, b_loc);
    f.save_to_file("OUT/filter.filt");

    Filter g(n, b_loc, 10);
    g.dump("OUT/bsft.filt");
    return 0;
}
