# polynomial/Makefile — build & gold regression for Khovanov Poincaré path
#
# Usage (from repository root):
#   make -C polynomial gold
#   make -C polynomial check
#   make -C polynomial clean
#
# Or from this directory after placing files:
#   make gold
#
# Produces:
#   build/test_gold_khovanov   — step-5 gold table
#   build/libkhpoly.a          — static lib (optional link target)

CC      ?= gcc
CFLAGS  ?= -std=c11 -Wall -Wextra -O2
LDFLAGS ?=

# Paths relative to this Makefile (lives in polynomial/)
ROOT    := ..
BUILD   := build

INCLUDES := \
	-I$(ROOT) \
	-I. \
	-Ihomology \
	-Ipolynomial \
	-Iencoder

# Production polynomial stack (no topology / CobMatrix required)
POLY_SRCS := \
	polynomial/BivariatePoly.c \
	homology/KhPoincare.c \
	homology/TopologyBridge.c \
	encoder/KnotEncoder.c \
	$(ROOT)/IntegerMatrix.c

POLY_OBJS := $(patsubst %.c,$(BUILD)/%.o,$(notdir $(POLY_SRCS)))

# Map object names (flat under build/)
$(BUILD)/BivariatePoly.o:     polynomial/BivariatePoly.c
$(BUILD)/KhPoincare.o:        homology/KhPoincare.c
$(BUILD)/TopologyBridge.o:    homology/TopologyBridge.c
$(BUILD)/KnotEncoder.o:       encoder/KnotEncoder.c
$(BUILD)/IntegerMatrix.o:     $(ROOT)/IntegerMatrix.c

TEST_GOLD_SRC := test/test_gold_khovanov.c
TEST_GOLD_BIN := $(BUILD)/test_gold_khovanov

.PHONY: all gold check lib clean dirs

all: gold

dirs:
	@mkdir -p $(BUILD)

# Pattern: compile any listed source into build/<basename>.o
$(BUILD)/BivariatePoly.o $(BUILD)/KhPoincare.o $(BUILD)/TopologyBridge.o \
$(BUILD)/KnotEncoder.o $(BUILD)/IntegerMatrix.o: | dirs
	$(CC) $(CFLAGS) $(INCLUDES) -c $(filter %.c,$^) -o $@

lib: $(BUILD)/libkhpoly.a

$(BUILD)/libkhpoly.a: $(BUILD)/BivariatePoly.o $(BUILD)/KhPoincare.o \
		$(BUILD)/TopologyBridge.o $(BUILD)/KnotEncoder.o $(BUILD)/IntegerMatrix.o
	ar rcs $@ $^

$(TEST_GOLD_BIN): $(TEST_GOLD_SRC) $(BUILD)/BivariatePoly.o $(BUILD)/KhPoincare.o \
		$(BUILD)/TopologyBridge.o $(BUILD)/KnotEncoder.o $(BUILD)/IntegerMatrix.o | dirs
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $(TEST_GOLD_SRC) \
		$(BUILD)/BivariatePoly.o $(BUILD)/KhPoincare.o \
		$(BUILD)/TopologyBridge.o $(BUILD)/KnotEncoder.o \
		$(BUILD)/IntegerMatrix.o $(LDFLAGS)

gold: $(TEST_GOLD_BIN)
	@echo "=== running gold tests ==="
	./$(TEST_GOLD_BIN)

check: gold
	@echo "=== ALL CHECKS PASSED ==="

clean:
	rm -rf $(BUILD)
