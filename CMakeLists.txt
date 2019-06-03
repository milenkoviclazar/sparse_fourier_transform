cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/include -L/usr/local/lib -Wall -O2 -Wconversion -std=c++03")

set(SRCS_OPT
        block_sparse/bsft.cc
        block_sparse/fftw.cc
        block_sparse/computefourier.cc
        block_sparse/utils.cc
        block_sparse/timer.cc
        block_sparse/filters.cc
        block_sparse/parameters.cc)

set(SRCS_MIT
        sFFT/computefourier.cc
        sFFT/fftw.cc
        sFFT/filters.cc
        sFFT/parameters.cc
        sFFT/plot.cc
        sFFT/timer.cc
        sFFT/utils.cc)

add_executable(bsft_tests block_sparse/bsft_tests.cc ${SRCS_OPT})
target_link_libraries(bsft_tests fftw3  m)

add_executable(precompute_filters block_sparse/precompute_filters.cc ${SRCS_OPT})
target_link_libraries(precompute_filters fftw3  m)

add_executable(generate_signals block_sparse/generate_signals.cc ${SRCS_OPT})
target_link_libraries(generate_signals fftw3  m)

add_executable(k_sparse_tuning block_sparse/k_sparse_tuning.cc ${SRCS_OPT})
target_link_libraries(k_sparse_tuning fftw3  m)

add_executable(run_experiment block_sparse/run_experiment.cc ${SRCS_OPT})
target_link_libraries(run_experiment fftw3  m)

add_executable(kaiserbessel sFFT/kaiserbessel.cc ${SRCS_MIT})
target_link_libraries(kaiserbessel fftw3 m)
