#ifndef KH_MEMORY_H
#define KH_MEMORY_H
/* Standalone solver policy: allocation failure terminates the computation.
 * Zero-size requests receive valid storage; no failed allocation is algebraic zero.
 * Define KH_TEST_ALLOC_FAILURE only for allocator-failure regression builds. */
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
static inline _Noreturn void kh_fatal(const char *message) {
    fprintf(stderr, "Error: %s\n", message);
    exit(EXIT_FAILURE);
}
static inline size_t kh_size_mul(size_t a, size_t b) {
    if (b && a > SIZE_MAX / b) kh_fatal("allocation size overflow");
    return a * b;
}
static inline void kh_allocation_check(void) {
#ifdef KH_TEST_ALLOC_FAILURE
    /* Per-translation-unit countdown, deliberately absent in release builds. */
    static long remaining = -2;
    if (remaining == -2) {
        const char *s = getenv("KH_FAIL_ALLOC_AFTER");
        remaining = s ? strtol(s, NULL, 10) : -1;
    }
    if (remaining == 0) kh_fatal("injected allocation failure");
    if (remaining > 0) --remaining;
#endif
}
static inline void *kh_malloc(size_t size) {
    kh_allocation_check();
    void *p = malloc(size ? size : 1);
    if (!p) kh_fatal("out of memory");
    return p;
}
static inline void *kh_calloc(size_t count, size_t size) {
    size_t bytes = kh_size_mul(count, size);
    kh_allocation_check();
    void *p = calloc(bytes ? bytes : 1, 1);
    if (!p) kh_fatal("out of memory");
    return p;
}
static inline void *kh_realloc(void *old, size_t size) {
    kh_allocation_check();
    void *p = realloc(old, size ? size : 1);
    if (!p) kh_fatal("out of memory");
    return p;
}
static inline int64_t kh_add64(int64_t a, int64_t b) {
    if ((b > 0 && a > INT64_MAX - b) || (b < 0 && a < INT64_MIN - b))
        kh_fatal("int64_t addition overflow");
    return a + b;
}
static inline int64_t kh_mul64(int64_t a, int64_t b) {
    if (a && b) {
        if ((a > 0 && b > 0 && a > INT64_MAX / b) ||
            (a > 0 && b < 0 && b < INT64_MIN / a) ||
            (a < 0 && b > 0 && a < INT64_MIN / b) ||
            (a < 0 && b < 0 && b < INT64_MAX / a))
            kh_fatal("int64_t multiplication overflow");
    }
    return a * b;
}
/* Reinterpret a wrapping 32-bit hash without an out-of-range signed cast. */
static inline int kh_hash_int(uint32_t bits) {
    int64_t value = bits <= INT32_MAX ? (int64_t)bits
                                    : (int64_t)bits - INT64_C(4294967296);
    return (int)value;
}
static inline int kh_int(int64_t value) {
    if (value < INT_MIN || value > INT_MAX) kh_fatal("integer size or grading overflow");
    return (int)value;
}
#endif
