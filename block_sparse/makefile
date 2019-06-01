RM=rm -f
CXXFLAGS=-O2 -std=c++03 -I ${FFTW_ROOT}/include -L ${FFTW_ROOT}/lib
LDLIBS=-lfftw3 -lfftw3f -lm

DESTS=precompute_filters bsft_tests generate_signals k_sparse_tuning run_experiment
SRCS:=$(filter-out $(DESTS:%=%.cc), $(wildcard *.cc))
OBJS:=$(subst .cc,.o,$(SRCS))

all: $(DESTS)

$(DESTS): $(OBJS)

clean:
	$(RM) *.o $(DESTS)
