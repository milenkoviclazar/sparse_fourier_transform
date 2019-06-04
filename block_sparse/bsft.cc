#include <iostream>
#include <map>
#include <cassert>
#include <cstring>
#include <fstream>
#include <set>
#include <algorithm>
#include <complex.h>
#include <iomanip>
#include <cmath>

#include "bsft.h"
#include "computefourier.h"
#include "parameters.h"
#include "utils.h"
#include "filters.h"

double TOTAL_TIME = 0.0;
clock_t prev;
#define start_timer() {prev = clock();}
#define stop_timer() {TOTAL_TIME += clock() - prev;}
#define get_curr_time(s) {cout << "\033[33mTIMER " << (s) << ": " << (TOTAL_TIME + (clock() - prev)) / CLOCKS_PER_SEC << "\033[0m\n";}
#define get_sample_count(s, n) {cout << "\033[34mSAMPLES " << (s) << ": " << 1.0 * count_samples_taken((n)) / n << "\033[0m\n";}

//#define RELEASE 1

bool *idx_taken;
complex_t **cache;
map<pair<int, int>, FilterHikp> filters;
map<pair<int, int>, FilterHikp> filters_est;
int *sbB_loc, *sbB_thresh, *sbB_est, *sbloc_loops, *sbest_loops, *sbthresh_loops;

using namespace std;

// ********************************************* FILTER ********************************************* //
void Filter::normalize(complex_t a[], int n) {
    complex_t m = 0;
    for (int i = 0; i < n; i++) {
        if (cabs(a[i]) > cabs(m)) {
            m = a[i];
        }
    }
    for (int i = 0; i < n; i++) {
        a[i] /= m;
    }
}

void Filter::convolve(complex_t a[], complex_t b[], int n) {
    for (int i = 0; i < n; i++) {
        x[i] = a[i];
        y[i] = b[i];
    }
    fftw_dft(x, n, x);
    fftw_dft(y, n, y);
    for (int i = 0; i < n; i++) {
        tmp[i] = x[i] * y[i];
    }
    fftw_dft(tmp, n, tmp, 1);
    normalize(tmp, n);
    for (int i = 0; i < n; i++) {
        a[i] = tmp[i];
    }
}

void Filter::resize(double energy_percent) {
    double total_energy = 0.0;
    for (int i = 0; i < n / 2; i++) {
        total_energy += cabs(time[i]);
    }
    double curr_energy = 0.0;
    for (size = 0; size < n / 2; size++) {
        curr_energy += cabs(time[size]);
        if (curr_energy >= energy_percent * total_energy) {
            break;
        }
    }
}

void Filter::dump(const char *filename) {
#ifdef RELEASE
    assert(1 == 0);
#endif
    ofstream ofs;
    ofs.open(filename);
    ofs << n << endl;
    for (int i = 0; i < n; i++) {
        ofs << creal(time[i]);
        if (i + 1 != n) {
            ofs << ",";
        }
    }
    ofs << endl;
    for (int i = 0; i < n; i++) {
        ofs << cimag(time[i]);
        if (i + 1 != n) {
            ofs << ",";
        }
    }
    ofs << endl;
    for (int i = 0; i < n; i++) {
        ofs << creal(freq[i]);
        if (i + 1 != n) {
            ofs << ",";
        }
    }
    ofs << endl;
    for (int i = 0; i < n; i++) {
        ofs << cimag(freq[i]);
        if (i + 1 != n) {
            ofs << ",";
        }
    }
    ofs << endl;
    ofs.close();
}

void Filter::dump_binary(const char *path) {
    char filename[100];
    sprintf(filename, "%s/%d_%d_%d.filt", path, n, B, F);
    cout << filename << endl;
    ofstream ofs(filename, ofstream::binary);
    int size = n * sizeof(complex_t);
    char *buff = new char[size];
    memcpy(buff, time, size);
    ofs.write(buff, size);
    memcpy(buff, freq, size);
    ofs.write(buff, size);
    ofs.close();
    delete[] buff;
}

Filter::Filter(const char *path, int n, int B, int F) {
    this->n = n;
    this->B = B;
    this->F = F;
    x = y = tmp = W = W_ = rect = 0;

    char filename[100];
    sprintf(filename, "%s/%d_%d_%d.filt", path, n, B, F);
    ifstream ifs(filename, ifstream::binary);
    if (ifs.bad() || ifs.fail()) {
        cerr << "input filter doesn't exist for " << n << ", " << B << ", " << F << endl;
        this->n = -1;
        this->B = -1;
        this->F = -1;
        return;
    }
    int size = n * sizeof(complex_t);
    char *buff = new char[size];
    ifs.read(buff, size);
    time = new complex_t[size];
    memcpy(time, buff, size);
    ifs.read(buff, size);
    freq = new complex_t[size];
    memcpy(freq, buff, size);

    ifs.close();
    delete[] buff;
}

complex_t * make_dolphchebyshev_t_edited(double lobefrac, double tolerance, int &w){
    w = min(w, int((1 / M_PI) * (1/lobefrac) * acosh(1./tolerance)));
    if (!(w%2))
        w--;
    complex_t *x = (complex_t *)malloc(w*sizeof(*x));
    double t0 = cosh(acosh(1/tolerance) / (w-1));
    for(int i = 0; i < w; i++){
        x[i] = Cheb(w-1, t0 * cos(M_PI * i / w)) * tolerance;
    }
    fftw_dft(x, w, x);
    shift(x, w, w/2);
    for(int i = 0; i < w; i++)
        x[i] = creal(x[i]);
    return x;
}

Filter::Filter(int n, int B, int F, bool dolphchebyshev) {
    this->n = n;
    this->B = B;
    this->F = F;
    if (dolphchebyshev) {
        double frac = 0.5 / B;
        int b = int(1.4 * 1.1 * 2 * frac * n); // TODO: taken from HIKP... shall be tuned
        complex_t *filtert3 = (complex_t *)calloc(n, sizeof(complex_t));
        complex_t *filterf3 = (complex_t *)calloc(n, sizeof(complex_t));

        int w3 = n;
        tmp = make_dolphchebyshev_t_edited(frac, 1e-9, w3);
        memcpy(filtert3, tmp, w3*sizeof(*filtert3)); free(tmp);
        make_multiple_t(filtert3, w3, n, b);
        for (int i = 0; i < n; i++) {
            filtert3[i] *= n;
        }
        shift(filtert3, n, - w3 / 2);
        fftw_dft(filterf3, n, filtert3); // from time to freq shall need a normalization ?
        for (int i = 0; i < n; i++) {
            filterf3[i] /= n;
        }
        time = filtert3;
        freq = filterf3;
        size = w3 / 4;
    } else {
        x = (complex_t *) calloc(n, sizeof(complex_t));
        y = (complex_t *) calloc(n, sizeof(complex_t));
        tmp = (complex_t *) calloc(n, sizeof(complex_t));
        W = (complex_t *) calloc(n, sizeof(complex_t));
        W_ = (complex_t *) calloc(n, sizeof(complex_t));
        rect = (complex_t *) calloc(n, sizeof(complex_t));
        time = (complex_t *) calloc(n, sizeof(complex_t));
        freq = (complex_t *) calloc(n, sizeof(complex_t));

        int Bprim = min(n >> 1, B << 4);
        for (int i = 0; i < Bprim / 2; i++) {
            W[i] = 1;
            W[(n - i) % n] = 1;
        }

        for (int i = 0; i < Bprim / 2; i++) {
            rect[i] = 1;
            rect[(n - i) % n] = 1;
        }
        for (int i = 0; i < F - 1; i++) {
            convolve(W, rect, n);
        }

        fftw_dft(W_, n, W);
        normalize(W_, n);

        for (int i = 0; i < n; i++) {
            rect[i] = 0;
        }
        for (int i = 0; i < 0.95 * n / B; i++) {
            rect[i] = 1;
            rect[(n - i) % n] = 1;
        }

        convolve(W_, rect, n);
        for (int i = 0; i < n; i++) {
            freq[i] = W_[i];
        }
        size = min(F * Bprim / 2, n / 2);
        fftw_dft(time, n, freq, 1); // from freq to time no normalization
    }
}

void Filter::clr() {
    if (x != 0) {
        free(x);
        free(y);
        free(tmp);
        free(W);
        free(W_);
        free(rect);
        x = y = tmp = W = W_ = rect = 0;
    }

    free(time);
    free(freq);
}


// ********************************************* UTILITIES ********************************************* //
double cabs2(complex_t x) {
    return (__real__ x) * (__real__ x) + (__imag__ x) * (__imag__ x);
}

int get_block_idx(int i, int k1, int n) {
    return ((i + k1 / 2 - 1) / k1) % (n / k1);
}

void get_interval_bounds(int j, int k1, int &l, int &r, const int n) {
    // TODO: this has a bug when the returned interval is (- k1/2, k1/2]
    // TODO: in that case l = 2^n - k1/2 and r = k1/2 !!!
    l = j - k1 / 2 + 1;
    while (l < 0) {
        l += n;
    }
    r = j + k1 / 2;
    while (r >= n) {
        r -= n;
    }
}

int round_to_power2(int x) {
    for (int i = 1; i <= 31; i++) {
        if (x <= (1 << i)) {
            return (1 << i);
        }
    }
    assert(false);
    return 0;
}

void save_signal(const char *filename, complex_t X[], int n) {
#ifdef RELEASE
    assert(1 == 0);
#endif
    ofstream ofs;
    ofs.open(filename);
    ofs << n << endl;
    for (int i = 0; i < n; i++) {
        ofs << __real__ X[i];
        if (i + 1 < n) {
            ofs << ", ";
        }
    }
    ofs << endl;
    for (int i = 0; i < n; i++) {
        ofs << __imag__ X[i];
        if (i + 1 < n) {
            ofs << ", ";
        }
    }
    ofs << endl;
    ofs.close();
}

void save_signal_binary(const char *path, complex_t X[], set<int> &blocks, int n, int k0, int k1, int idx) {
    char filename[256];
    sprintf(filename, "%srandom_%d_%d_%d__%d.in", path, n, k0, k1, idx);
    ofstream ofs(filename, ofstream::binary);
    size_t size = n * sizeof(complex_t);
    char *buff = new char[size];
    memcpy(buff, X, size);
    ofs.write(buff, size);
    int *tmp = new int[k0];
    int _idx = 0;
    for (set<int>::iterator it = blocks.begin(); it != blocks.end(); it++) {
        tmp[_idx++] = *it;
    }
    size = k0 * sizeof(int);
    memcpy(buff, tmp, size);
    ofs.write(buff, size);

    delete[] buff;
    ofs.close();
}

bool load_signal_binary(const char *path, complex_t X[], set<int> &blocks, int n, int k0, int k1, int idx) {
    char filename[256];
    sprintf(filename, "%s/random_%d_%d_%d__%d.in", path, n, k0, k1, idx);
    ifstream ifs(filename, ifstream::binary);
    if (ifs.bad() || ifs.fail()) {
        cerr << "input signal doesn't exist for " << n << ", " << k0 << ", " << k1 << ", " << idx << endl;
        return false;
    }

    size_t size = n * sizeof(complex_t);
    char *buff = new char[size];
    ifs.read(buff, size);
    memcpy(X, buff, size);
    size = k0 * sizeof(int);
    ifs.read(buff, size);
    int *tmp = new int[k0];
    memcpy(tmp, buff, size);
    for (int i = 0; i < k0; i++) {
        blocks.insert(tmp[i]);
    }
    delete[] buff;
    ifs.close();
    return true;
}

void initialize_globals(int n, int k1) {
    idx_taken = new bool[n];
    memset(idx_taken, false, n * sizeof(bool));
    cache = new complex_t *[2 * k1];
    for (int i = 0; i < 2 * k1; i++) {
        cache[i] = new complex_t[n / k1];
        memset(cache[i], 0xfe, n / k1 * sizeof(complex_t));
    }
}

void destruct_globals(int n, int k1) {
    delete[] idx_taken;
    for (int i = 0; i < 2 * k1; i++) {
        delete[] cache[i];
    }
    delete[] cache;
}

int count_samples_taken(int n) {
    int samples_taken = 0;
    for (int i = 0; i < n; i++) {
        if (idx_taken[i]) {
            samples_taken++;
        }
    }
    return samples_taken;
}

// ********************************************* BLOCK LOCATION ********************************************* //
complex_t get_downsampled(complex_t *X, Filter &G, int r, int j, int k1, int n) {
    if (*(unsigned long long *) (&cache[r][j]) != 0xfefefefefefefefe) {
        return cache[r][j];
    }
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
            result += G.time[idxG] * take_sample(X, idxX);
        }
    }
    cache[r][j] = result / k1;
    return cache[r][j];
}

void hash_reduced(
        int n,
        int k1,
        int r,
        complex_t *out,
        complex_t *X,
        long long B,
        Filter &G_hash,
        Filter &G_downsample,
        long long sigma,
        long long delta) {
    int nk1 = n / k1;
    long long count = G_hash.size / B;
    for (long long j = 0; j < B; j++) {
        long long Bi = nk1 - j - count * B / 2;
        while (Bi < 0) {
            Bi += nk1;
        }
        Bi = (Bi / B) * B;
        out[j] = 0;
        for (int iter = 0; iter < count; iter++) {
            long long idx = sigma * (delta + j + Bi) % nk1;
            out[j] += get_downsampled(X, G_downsample, r, idx, k1, n) * G_hash.time[(j + Bi) % nk1];
            Bi += B;
            while (Bi >= nk1) {
                Bi -= nk1;
            }
        }
        out[j] *= (double) B / (double) nk1;
    }
}

double hash_to_bins_reduced(
        int n,
        int k1,
        int r,
        complex_t *X,
        complex_t *U,
        complex_t *U_,
        Filter &G_hash,
        Filter &G_downsample,
        int B,
        int sigma,
        int delta) {
    hash_reduced(n, k1, r, U, X, B, G_hash, G_downsample, sigma, delta);
    fftw_dft(U_, B, U);
    double ret = 0.0;
    for (int i = 0; i < B; i++) {
        ret += cabs2(U_[i] / B);
    }
    return ret;
}

map<int, complex_t> sparse_fft_wrapper(
        complex_t *X,
        int n,
        int k,
        int sB_est,
        int sB_thresh,
        int sB_loc,
        int sW_Comb,
        int sComb_loops,
        int sthreshold_loops,
        int sloc_loops,
        int sest_loops) {
    return outer_loop(X, n, filters[make_pair(n, k)], filters_est[make_pair(n, k)], sB_est, sB_thresh, sB_loc, sW_Comb,
                      sComb_loops, sthreshold_loops, sloc_loops, sloc_loops + sest_loops);
}

void set_sparse_fft_parameters(
        int n,
        int k,
        int &B_est,
        int &B_loc,
        int &B_thresh,
        int &W_Comb,
        int &loc_loops,
        int &est_loops,
        int &threshold_loops,
        int &Comb_loops,
        bool recreateFilters,
        double Bcst_loc,
        double Bcst_est,
        bool getDefaultValues) {
    double Comb_cst = 16;
    double tolerance_loc = 0.01;
    double tolerance_est = 0.01;

    if (getDefaultValues) {
        Bcst_loc = 1.0;
        Bcst_est = 1.0;
        if (k <= 2) {
            loc_loops = 6;
            est_loops = 8;
        } else if (k <= 4) {
            loc_loops = 7;
            est_loops = 10;
        } else if (k <= 8) {
            loc_loops = 8;
            est_loops = 10;
        } else if (k <= 16) {
            loc_loops = 8;
            est_loops = 12;
        } else {
            loc_loops = 9;
            est_loops = 12;
        }
        threshold_loops = loc_loops - 1;
    }


    real_t BB_loc = (unsigned) (Bcst_loc * sqrt((double) n * k / (log2(n))));
    real_t BB_est = (unsigned) (Bcst_est * sqrt((double) n * k / (log2(n))));

    double lobefrac_loc = 0.5 / (BB_loc);
    double lobefrac_est = 0.5 / (BB_est);

    int b_loc = int(1.2 * 1.1 * ((double) n / BB_loc));
    int b_est = int(1.4 * 1.1 * ((double) n / BB_est));

    B_loc = floor_to_pow2(BB_loc);
    B_thresh = 2 * k;
    B_est = floor_to_pow2(BB_est);

    W_Comb = floor_to_pow2(Comb_cst * n / B_loc);

    if (recreateFilters || filters.find(make_pair(n, k)) == filters.end()) {
        int w_loc;
        complex_t *filtert = make_dolphchebyshev_t(lobefrac_loc, tolerance_loc, w_loc);
        filters[make_pair(n, k)] = make_multiple_t(filtert, w_loc, n, b_loc);
    }

    if (recreateFilters || filters_est.find(make_pair(n, k)) == filters_est.end()) {
        int w_est;
        complex_t *filtert_est = make_dolphchebyshev_t(lobefrac_est, tolerance_est, w_est);
        filters_est[make_pair(n, k)] = make_multiple_t(filtert_est, w_est, n, b_est);
    }
}

bool cmp_abs(pair<int, complex_t> a, pair<int, complex_t> b) {
    return cabs2(a.second) > cabs2(b.second);
}

void multi_block_locate(
        complex_t X[],
        int n,
        int k0,
        int k1,
        Filter &H_hash,
        Filter &H_downsample,
        set<int> &blocks,
        int B,
        int iter_loc,
        int iter_budget) {
    int qlim = 0;
    while ((1 << qlim) <= k0) {
        qlim++;
    }
    qlim++;
    double *gamma = (double *) calloc(2 * k1, sizeof(double));
    int *s = (int *) calloc(2 * k1, sizeof(int));
    memset(s, -1, 2 * k1 * sizeof(int));
    double *gamma_cumulative = (double *) calloc(2 * k1, sizeof(double));
    double *q_cumulative = (double *) calloc(qlim, sizeof(double));
    complex_t *U = (complex_t *) calloc(B, sizeof(complex_t));
    complex_t *U_ = (complex_t *) calloc(B, sizeof(complex_t));

    double sum = 0.0;
    q_cumulative[0] = 0.5;
    for (int i = 1; i < qlim; i++) {
        q_cumulative[i] = q_cumulative[i - 1] + (1.0 / (1LL << (i + 1)));
        sum = q_cumulative[i];
    }
    for (int i = 0; i < qlim; i++) {
        q_cumulative[i] /= sum;
    }
    while (iter_loc--) {
        int Delta = rand_int(n / k1);
        int sigma = (rand_int(n / k1 / 2) * 2) + 1;
        sum = 0.0;
        for (int i = 0; i < 2 * k1; i++) {
            gamma[i] = hash_to_bins_reduced(n, k1, i, X, U, U_, H_hash, H_downsample, B, sigma, Delta);
            sum += gamma[i];
        }
        gamma_cumulative[0] = gamma[0] / sum;
        for (int i = 1; i < 2 * k1; i++) {
            gamma_cumulative[i] = gamma_cumulative[i - 1] + gamma[i] / sum;
        }
        for (int _iter = 0; _iter < iter_budget; _iter++) {
            double rnd = drand48();
            int q, r;
            for (r = 0; r < 2 * k1; r++) {
                if (gamma_cumulative[r] >= rnd) {
                    break;
                }
            }
            rnd = drand48();
            for (q = 0; q < qlim; q++) {
                if (q_cumulative[q] >= rnd) {
                    break;
                }
            }
            if (r == 2 * k1) {
                continue;
            }
            s[r] = max(s[r], q);
        }
    }
    for (int i = 0; i < 2 * k1; i++) {
        if (s[i] == -1) {
            continue;
        }
        int idx = s[i];
        int sparsity = (1 << idx);
        if (sparsity > k0) {
            sparsity = k0;
        }
        int unusedVariable = 0;
        map<int, complex_t> ans = outer_loop(X, n / k1, filters[make_pair(n / k1, sparsity)],
                                             filters_est[make_pair(n / k1, sparsity)], sbB_est[idx], sbB_thresh[idx],
                                             sbB_loc[idx], unusedVariable,
                                             unusedVariable, sbthresh_loops[idx], sbloc_loops[idx],
                                             sbloc_loops[idx] + sbest_loops[idx],
                                             H_downsample, i, k1);
        vector<pair<int, complex_t> > v;
        for (map<int, complex_t>::iterator it = ans.begin(); it != ans.end(); it++) {
            v.push_back(*it);
//            if (cabs(it->second) > 0.01) {
//                blocks.insert(it->first * k1);
//            }
        }
        sort(v.begin(), v.end(), cmp_abs);
        for (vector<pair<int, complex_t> >::iterator it = v.begin(); it != v.end(); it++) {
            if (cabs(it->second) > 0.1 || blocks.size() < k0) {
                blocks.insert(it->first * k1);
            }
        }
    }
    free(gamma);
    free(s);
    free(gamma_cumulative);
    free(q_cumulative);
    free(U_);
    free(U);
}

// ********************************************* VALUE ESTIMATION ********************************************* //
void hash_normal(
        int n,
        complex_t *out,
        complex_t *X,
        long long B,
        Filter &G,
        long long sigma,
        long long delta) {
    long long L = G.size;
    long long count = L / B;
    for (long long j = 0; j < B; j++) {
        long long Bi = n - j - L / 2;
        while (Bi < 0) {
            Bi += n;
        }
        Bi = (Bi / B) * B;
        out[j] = 0;
        for (int iter = 0; iter < count; iter++) {
            long long idx = sigma * (delta + j + Bi) % n;
            out[j] += take_sample(X, idx) * G.time[(j + Bi) % n];
            Bi += B;
            while (Bi >= n) {
                Bi -= n;
            }
        }
        out[j] *= (double) B / (double) n;
    }
}

void hash_to_bins(
        int n,
        complex_t *X,
        complex_t *U,
        complex_t *U_,
        Filter &G,
        int B,
        int sigma,
        int Delta) {
    hash_normal(n, U, X, B, G, sigma, Delta);
    fftw_dft(U_, B, U);
    for (int i = 0; i < B; i++) {
        U_[i] /= B;
    }
}

bool cmp_real(const complex_t &A, const complex_t &B) {
    return creal(A) < creal(B);
}

bool cmp_imag(const complex_t &A, const complex_t &B) {
    return cimag(A) < cimag(B);
}

long long h(long long j, long long B, long long n, long long sigma) {
    long long ret = sigma * j * B / n;
    return ret % B;
}

set<int> prune_location(
        complex_t X[],
        int n,
        int k0,
        int k1,
        const set<int> &blocks,
        Filter &G,
        int B,
        int T,
        double theta) {
    int *freqs = (int *) calloc(blocks.size() * k1, sizeof(complex_t));
    complex_t *U = new complex_t[B];
    complex_t *U_ = new complex_t[B];
    int freq_cnt = 0;
    for (set<int>::iterator it = blocks.begin(); it != blocks.end(); it++) {
        int lo, hi;
        get_interval_bounds(*it, k1, lo, hi, n);
        for (int i = lo; i <= hi; i++) {
            freqs[freq_cnt++] = i;
        }
    }
    double **power = new double *[blocks.size()];
    for (int i = 0; i < blocks.size(); i++) {
        power[i] = new double[T];
        memset(power[i], 0, T * sizeof(double));
    }
    const complex_t OMEGA_N = cexp(2 * I * M_PI / n);

    for (int t = 0; t < T; t++) {
        int Delta = rand_int(n / k1);
        int sigma = (rand_int(n / k1 / 2) * 2) + 1;
        hash_to_bins(n, X, U, U_, G, B, sigma, Delta);
        for (int i = 0; i < freq_cnt; i++) {
            int f = freqs[i];
            complex_t tmp = U_[h(f, B, n, sigma)] * cpow(OMEGA_N, -sigma * Delta * f);
            power[i / k1][t] += cabs2(tmp);
        }
    }
    set<int> ret;
    vector<pair<int, complex_t> > v;
    for (int i = 0; i < freq_cnt; i++) {
        int ik1 = i / k1;
        nth_element(power[ik1], power[ik1] + (T / 2), power[ik1] + T);
        v.push_back(make_pair(get_block_idx(freqs[i], k1, n) * k1, power[ik1][T / 2]));

//        if (power[ik1][T / 2] < theta) {
//            continue;
//        }
//        ret.insert(get_block_idx(freqs[i], k1, n) * k1);
    }
    sort(v.begin(), v.end(), cmp_abs);
//    cout << v.size() << " " << k0 <<endl;
    for (int i = 0; i < v.size(); i++) {
        if (ret.size() < k0) {
            ret.insert(v[i].first);
        }
    }
//    cout << ret.size() << " " << k0 << endl;
    delete[] U;
    delete[] U_;
    for (int i = 0; i < blocks.size(); i++) {
        delete[] power[i];
    }
    delete[] power;
    return ret;
};

void estimate_values(
        complex_t X[],
        int n,
        int k0,
        int k1,
        const set<int> &blocks,
        map<int, complex_t> &result,
        Filter &G,
        int B,
        int T,
        double theta) {
    int *freqs = (int *) calloc(blocks.size() * k1, sizeof(complex_t));
    complex_t *U = new complex_t[B];
    complex_t *U_ = new complex_t[B];
    int freq_cnt = 0;
    for (set<int>::iterator it = blocks.begin(); it != blocks.end(); it++) {
        int lo, hi;
        get_interval_bounds(*it, k1, lo, hi, n);
        for (int i = lo; i <= hi; i++) {
            freqs[freq_cnt++] = i;
        }
    }
    complex_t **W = new complex_t *[freq_cnt];
    for (int i = 0; i < freq_cnt; i++) {
        W[i] = new complex_t[T];
        memset(W[i], 0, T * sizeof(complex_t));
    }
    double **power = new double *[blocks.size()];
    for (int i = 0; i < blocks.size(); i++) {
        power[i] = new double[T];
        memset(power[i], 0, T * sizeof(double));
    }
    const complex_t OMEGA_N = cexp(2 * I * M_PI / n);

    for (int t = 0; t < T; t++) {
        int Delta = rand_int(n / k1);
        int sigma = (rand_int(n / k1 / 2) * 2) + 1;
        hash_to_bins(n, X, U, U_, G, B, sigma, Delta);
        for (int i = 0; i < freq_cnt; i++) {
            int f = freqs[i];
            complex_t tmp = U_[h(f, B, n, sigma)] * cpow(OMEGA_N, -sigma * Delta * f);
            W[i][t] = tmp;
            power[i / k1][t] += cabs2(tmp);
        }
    }
    for (int i = 0; i < freq_cnt; i++) {
        int ik1 = i / k1;
        nth_element(power[ik1], power[ik1] + (T / 2), power[ik1] + T);
        if (power[ik1][T / 2] < theta) {
            continue;
        }
        nth_element(W[i], W[i] + (T / 2), W[i] + T, cmp_real);
        complex_t C;
        __real__ C = creal(W[i][T / 2]);
        nth_element(W[i], W[i] + (T / 2), W[i] + T, cmp_imag);
        __imag__ C = cimag(W[i][T / 2]);
        result[freqs[i]] = C;
    }
    delete[] U;
    delete[] U_;
    for (int i = 0; i < freq_cnt; i++) {
        delete[] W[i];
    }
    delete[] W;
    for (int i = 0; i < blocks.size(); i++) {
        delete[] power[i];
    }
    delete[] power;
};

// ********************************************* INTERFACE ********************************************* //
map<int, complex_t> bsft(
        complex_t X[],
        int n,
        int k0,
        int k1,
        Filter &G_val,
        Filter &H_downsample,
        Filter &H_hash,
        int B_loc,
        int B_val,
        int iter_loc,
        int iter_budget,
        int iter_val,
        double theta_val) {
    set<int> blocks;
    map<int, complex_t> result;
    multi_block_locate(X, n, k0, k1, H_hash, H_downsample, blocks, B_loc, iter_loc, iter_budget);
//    return prune_location(X, n, k0, k1, blocks, G_val, B_val, iter_val, theta_val);
    estimate_values(X, n, k0, k1, blocks, result, G_val, B_val, iter_val, theta_val);
//    return result;
}

set<int> bsft_location(
        complex_t X[],
        int n,
        int k0,
        int k1,
        Filter &G_val,
        Filter &H_downsample,
        Filter &H_hash,
        int B_loc,
        int B_val,
        int iter_loc,
        int iter_budget,
        int iter_val,
        double theta_val) {
    set<int> blocks;
    multi_block_locate(X, n, k0, k1, H_hash, H_downsample, blocks, B_loc, iter_loc, iter_budget);
    if (iter_val) {
        return prune_location(X, n, k0, k1, blocks, G_val, B_val, iter_val, theta_val);
    } else {
        return blocks;
    }

}

void set_parameters(
        const int n,
        const int k0,
        const int k1) {
    int qlim = 0;
    while ((1 << qlim) <= k0) {
        qlim++;
    }
    qlim++;
    if (sbB_loc != 0) {
        delete[] sbB_loc;
        delete[] sbB_thresh;
        delete[] sbB_est;
        delete[] sbloc_loops;
        delete[] sbest_loops;
        delete[] sbthresh_loops;
        sbB_loc = 0;
        sbB_thresh = 0;
        sbB_est = 0;
        sbloc_loops = 0;
        sbest_loops = 0;
        sbthresh_loops = 0;
    }
    sbB_loc = new int[qlim];
    sbB_thresh = new int[qlim];
    sbB_est = new int[qlim];
    sbloc_loops = new int[qlim];
    sbest_loops = new int[qlim];
    sbthresh_loops = new int[qlim];
    for (int i = 0; i < qlim; i++) {
        int sparsity = 1 << i;
        if (sparsity > k0) {
            sparsity = k0;
        }
        int unusedVar;

        set_sparse_fft_parameters(n / k1, sparsity, sbB_est[i], sbB_loc[i], sbB_thresh[i], unusedVar, sbloc_loops[i],
                                  sbest_loops[i], sbthresh_loops[i], unusedVar, false, 1.0, 1.0, true);
    }
}
