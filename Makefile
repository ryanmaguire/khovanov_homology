# Build and run the Khovanov homology scanner.
#
# Requirements:
#   Build: a C11 compiler and make.
#   Batch: zsh and the revised run_batch.sh beside this Makefile.
#
# Commands:
#
#   make : Build FullScanning; rebuild when sources or headers change.
#
#   make batch INPUT=knots.txt OUTPUT=results.csv JOBS=4
#   (Build if needed, then process DT codes with four workers.)
#
#   make batch INPUT=knots.txt OUTPUT=results.csv FORCE=1
#   (Explicitly replace an existing output CSV.)
#
# Input: one alphabetical DT code per line.
# Blank lines and # comments are allowed.
#
# Output CSV columns: dt,polynomial,exit_code.
# The polynomial records the integral free part.
# The batch runner also saves per-entry stderr logs.

CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

CORE = topology/Cap.c topology/CannedCobordismImpl.c topology/LCCC.c \
       topology/CobMatrix.c topology/Komplex.c topology/TangleKomplex.c \
       topology/TangleCompose.c topology/PDScanner.c topology/DTtoPD.c \
       IntegerMatrix.c polynomial/BivariatePoly.c

HEADERS = $(wildcard topology/*.h polynomial/*.h) IntegerMatrix.h KhMemory.h

INPUT ?=
OUTPUT ?= khovanov_batch.csv
JOBS ?= 1
FORCE ?= 0
export INPUT OUTPUT JOBS FORCE

.PHONY: all batch help clean
.DELETE_ON_ERROR:

# Default: build the scanner without running it.
all: FullScanning

FullScanning: topology/FullScanning.c $(CORE) $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) topology/FullScanning.c $(CORE) \
		$(LDFLAGS) $(LDLIBS) -o $@

# The runner handles file input, parallel workers, CSV output, and logs.
# Existing output is protected unless FORCE=1 is supplied.
batch: FullScanning
	@test -n "$$INPUT" || { \
		echo 'Usage: make batch INPUT=knots.txt [OUTPUT=results.csv] [JOBS=4] [FORCE=1]' >&2; \
		exit 2; \
	}
	@command -v zsh >/dev/null 2>&1 || { \
		echo 'Error: zsh is required by run_batch.sh.' >&2; \
		exit 2; \
	}
	@test -r ./run_batch.sh || { \
		echo 'Error: run_batch.sh must be beside this Makefile.' >&2; \
		exit 2; \
	}
	@set -- -f "$$INPUT" -b ./FullScanning -o "$$OUTPUT" -j "$$JOBS"; \
	case "$$FORCE" in \
		0) ;; \
		1) set -- "$$@" -F ;; \
		*) echo 'Error: FORCE must be 0 or 1.' >&2; exit 2 ;; \
	esac; \
	zsh ./run_batch.sh "$$@"

help:
	@printf '%s\n' \
		'make                       Build FullScanning.' \
		'make batch INPUT=knots.txt  Build if needed, then run a batch.' \
		'  OUTPUT=results.csv        Output path (default: khovanov_batch.csv).' \
		'  JOBS=4                    Parallel workers (default: 1).' \
		'  FORCE=1                   Allow replacement of an existing CSV.' \
		'make clean                 Remove the FullScanning executable.'

# Preserve input files, CSV results, and batch stderr logs.
clean:
	rm -f FullScanning