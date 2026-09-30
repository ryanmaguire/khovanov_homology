# polynomial/Makefile — matches FLAT layout (all sources in polynomial/)
#
# From repo root:
#   make -C polynomial gold
#   make -C polynomial low
#   make -C polynomial check
#   make -C polynomial clean

CC      ?= gcc
CFLAGS  ?= -std=c11 -Wall -Wextra -O2
LDFLAGS ?=

ROOT    := ..
BUILD   := build

INCLUDES := -I$(ROOT) -I.

POLY_OBJS := \
	$(BUILD)/BivariatePoly.o \
	$(BUILD)/KhPoincare.o \
	$(BUILD)/TopologyBridge.o \
	$(BUILD)/KnotEncoder.o \
	$(BUILD)/IntegerMatrix.o

.PHONY: all gold low check lib clean dirs

all: gold

dirs:
	@mkdir -p $(BUILD)

$(BUILD)/BivariatePoly.o: BivariatePoly.c BivariatePoly.h | dirs
	$(CC) $(CFLAGS) $(INCLUDES) -c BivariatePoly.c -o $@

$(BUILD)/KhPoincare.o: KhPoincare.c KhPoincare.h KnotEncoder.h BivariatePoly.h | dirs
	$(CC) $(CFLAGS) $(INCLUDES) -c KhPoincare.c -o $@

$(BUILD)/TopologyBridge.o: TopologyBridge.c TopologyBridge.h KhPoincare.h | dirs
	$(CC) $(CFLAGS) $(INCLUDES) -c TopologyBridge.c -o $@

$(BUILD)/KnotEncoder.o: KnotEncoder.c KnotEncoder.h BivariatePoly.h KhPoincare.h | dirs
	$(CC) $(CFLAGS) $(INCLUDES) -c KnotEncoder.c -o $@

$(BUILD)/IntegerMatrix.o: $(ROOT)/IntegerMatrix.c $(ROOT)/IntegerMatrix.h | dirs
	$(CC) $(CFLAGS) $(INCLUDES) -c $(ROOT)/IntegerMatrix.c -o $@

lib: $(BUILD)/libkhpoly.a

$(BUILD)/libkhpoly.a: $(POLY_OBJS)
	ar rcs $@ $^

$(BUILD)/test_gold_khovanov: test_gold_khovanov.c $(POLY_OBJS) | dirs
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ test_gold_khovanov.c $(POLY_OBJS) $(LDFLAGS)

$(BUILD)/test_low_crossing: test_low_crossing_javakh.c $(POLY_OBJS) | dirs
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ test_low_crossing_javakh.c $(POLY_OBJS) $(LDFLAGS)

gold: $(BUILD)/test_gold_khovanov
	@echo "=== running gold tests ==="
	./$(BUILD)/test_gold_khovanov

low: $(BUILD)/test_low_crossing
	@echo "=== running 4/5/6 crossing checks ==="
	./$(BUILD)/test_low_crossing

check: gold
	@echo "=== ALL CHECKS PASSED ==="

clean:
	rm -rf $(BUILD)
