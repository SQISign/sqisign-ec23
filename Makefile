CC=clang
CFLAGS=-Wall -Wextra -pedantic -std=gnu99 -I./include -I/usr/local/include
CFLAGSLINK=-lpari -lm -L/usr/local/lib/
DEBUG_FLAGS=-g
BENCH_FLAGS=-DNDEBUG -O3 -Os -march=native -mtune=native

default: all

# The implemented primes
PRIME_LIST=p6983 p3923
PRIME?=p3923

ASM=1  # Optimized assembly implementation used by default

ifndef ASM
EXTRA_MACROS=-DUSE_GMP
EXTRA_LIBS=-lgmp
endif

# The build directories
BUILDDIR = build/$(PRIME)
OBJDIR = $(BUILDDIR)/obj
TESTDIR = $(BUILDDIR)/test
BENCHDIR = $(BUILDDIR)/bench

$(BUILDDIR) $(OBJDIR) $(TESTDIR) $(BENCHDIR):
	mkdir -p $@

# Search paths
vpath %.h include
vpath %.c src:src/$(PRIME)
vpath %.s src:src/$(PRIME)
vpath %.o $(OBJDIR)

# Object files
$(OBJDIR)/%.o: %.c | $(OBJDIR)
	$(CC) $< $(CFLAGS) $(DEBUG_FLAGS) $(EXTRA_MACROS) -c -o $@

$(OBJDIR)/%.asm.o: %.s | $(OBJDIR)
	$(CC) $< -c -o $@


# Tests
$(TESTDIR)/%: test/%.c | $(TESTDIR)
	$(CC) $^ $(CFLAGS) $(DEBUG_FLAGS) $(CFLAGSLINK) $(EXTRA_LIBS) -o $@

# Benchmarks
$(BENCHDIR)/%: bench/%.c | $(BENCHDIR)
	$(CC) $^ $(CFLAGS) $(BENCH_FLAGS) $(EXTRA_MACROS) $(CFLAGSLINK) $(EXTRA_LIBS) -o $@

# Precomputations
$(BUILDDIR)/precomp $(BUILDDIR)/precomp_fp_const: | $(BUILDDIR)
	$(CC) $^ $(CFLAGS) $(DEBUG_FLAGS) $(CFLAGSLINK) $(EXTRA_LIBS) -o $@

precompute_fp: $(BUILDDIR)/precomp_fp_const
	$< > src/$(PRIME)/fp_const_generated.c

src/$(PRIME)/fp_const_generated.c: $(BUILDDIR)/precomp_fp_const
	$(MAKE) precompute_fp

precompute: $(BUILDDIR)/precomp
	$< > src/$(PRIME)/precomputed_generated.c

src/$(PRIME)/precomputed_generated.c: $(BUILDDIR)/precomp
	$(MAKE) precompute

# Velusqrt Tuning

$(BUILDDIR)/tunecycles: | $(BUILDDIR)
	$(CC) $^ $(CFLAGS) $(BENCH_FLAGS) $(EXTRA_MACROS) $(CFLAGSLINK) $(EXTRA_LIBS) -o $@

src/$(PRIME)/tunecycles.out: $(BUILDDIR)/tunecycles
	# 1 minute on 1.9GHz Kaby Lake
	time ./$< > $@

tune: src/tune2c src/$(PRIME)/tunecycles.out
	./src/tune2c < src/$(PRIME)/tunecycles.out > src/$(PRIME)/steps_tunecycles.c


# Run tests
TESTS=$(patsubst test/%.c,test_%,$(wildcard test/*.c))

$(TESTS:%=%): test_%: $(TESTDIR)/%
	@echo
	./$<

# Run benchmarks
BENCHS=$(patsubst bench/%.c,bench_%,$(wildcard bench/*.c))

$(BENCHS): bench_%: $(BENCHDIR)/%
	./$< >> $*.tsv

# Phony targets
check: $(TESTS)
benchmark: $(BENCHS)
tests: $(TESTS:test_%=$(TESTDIR)/%)
benchs: $(BENCHS:bench_%=$(BENCHDIR)/%)
all: tests benchs $(BUILDDIR)/tunecycles

distclean:
	rm -r build

.PHONY: distclean all tests benchs benchmark check tune default precompute


################################
## Object deps

ifdef ASM
ARITH_O		= uint.asm.o fp.asm.o
ARITH_SRC	= uint.s fp.s
else
ARITH_O		= uint.o fp.o
ARITH_SRC	= uint.c fp.c
endif
fp		= $(ARITH_O) gentobig.o rng.o fp2.o constants.o fp_const_generated.o
mont		= $(fp) mont.o poly.o curve.o tedwards.o steps.o steps_tunecycles.o
montxy		= $(mont)
isogenies	= $(mont) isogenies.o
isogenies_mult	= $(isogenies) toolbox.o
isom		= $(mont) isomorphism.o
two_walks	= $(isom) two_walks.o mitm.o
mitm		= $(two_walks)
mitm2 		= $(mitm) toolbox.o
bidim		= $(isogenies) $(two_walks)
quaternion_tools= $(mont) quaternion_tools.o ideal.o precomputed.o precomputed_generated.o toolbox.o
arith		= $(quaternion_tools)
klpt		= $(bidim) $(quaternion_tools) klpt.o
idiso		= $(klpt) idiso.o sqisign.o
sqisign         = $(idiso) verif.o

OBJECTS		= $(sqisign)

allobjects: $(addprefix $(OBJDIR)/, $(OBJECTS) )

$(BUILDDIR)/precomp_fp_const: precomp_fp_const.c $(ARITH_O) rng.o fp2.o constants.o
$(BUILDDIR)/precomp: precomp.c $(ARITH_O) gentobig.o rng.o fp2.o constants.o fp_const_generated.o mont.o curve.o tedwards.o poly.o steps.o steps_tunecycles.o
$(BUILDDIR)/tunecycles: tunecycles.c $(ARITH_SRC) isogenies.c mont.c fp2.c constants.c rng.c poly.c fp_const_generated.c curve.c tedwards.c steps.c steps_default.c

$(BENCHS:bench_%=$(BENCHDIR)/%): $(ARITH_SRC) $(filter-out %.asm.c,$(OBJECTS:%.o=%.c))

.SECONDEXPANSION:
$(TESTS:test_%=$(TESTDIR)/%): $(TESTDIR)/%: $$($$*)
