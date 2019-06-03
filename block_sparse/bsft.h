#ifndef BSFT_H
#define BSFT_H

#include <set>
#include <map>

#include "fftw.h"

typedef double _Complex complex_t;
#define take_sample(X, i) (idx_taken[i] = true, (X)[i])
#define rand_int(n) ((int)(drand48() * (n)))
extern bool *idx_taken;

using namespace std;

struct FilterHikp;

struct Filter {
    complex_t *time;
    complex_t *freq;
    int size;
    int n, B, F;
    complex_t *tmp;
    complex_t *x;
    complex_t *y;
    complex_t *W;
    complex_t *W_;
    complex_t *rect;

    void normalize(complex_t a[], int n);

    void convolve(complex_t a[], complex_t b[], int n);

    void resize(double energy_percent);

    void dump(const char *filename);

    void dump_binary(const char *path);

    Filter(int n, int B, int F, bool dolphchebyshev = true);

    Filter(const char *path, int n, int B, int F);

    void clr();

};

void get_interval_bounds(int j, int k1, int &l, int &r, int n);

void save_signal(const char *filename, complex_t X[], int n);

void set_parameters(
        int n,
        int k0,
        int k1);

void hash_to_bins(
        int n,
        complex_t *X,
        complex_t *U,
        complex_t *U_,
        Filter &G,
        int B,
        int sigma,
        int Delta);

int count_samples_taken(int n);

double cabs2(complex_t x);

int round_to_power2(int x);

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
        int delta);

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
        bool recreateFilters = true,
        double Bcst_loc = 1.0,
        double Bcst_est = 1.0,
        bool getDefaultValues = true);

complex_t get_downsampled(complex_t *X, Filter &G, int r, int j, int k1, int n);

void initialize_globals(int n, int k1);

void destruct_globals(int n, int k1);

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
        int sest_loops);

void save_signal_binary(const char *path, complex_t X[], set<int> &blocks, int n, int k0, int k1, int idx);
bool load_signal_binary(const char *path, complex_t X[], set<int> &blocks, int n, int k0, int k1, int idx);

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
        double theta_val);

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
        double theta_val);

set<int> prune_location(
        complex_t X[],
        int n,
        int k0,
        int k1,
        const set<int> &blocks,
        Filter &G,
        int B,
        int T,
        double theta);

int get_block_idx(int i, int k1, int n);

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
        double theta);

#endif