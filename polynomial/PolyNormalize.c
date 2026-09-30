#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct {
    long long coeff;
    long long q;
    long long t;
} Term;

static int cmp_terms(const void *a, const void *b) {
    const Term *ta = (const Term *)a;
    const Term *tb = (const Term *)b;
    if (ta->t != tb->t) return ta->t < tb->t ? -1 : 1;
    if (ta->q != tb->q) return ta->q < tb->q ? -1 : 1;
    if (ta->coeff != tb->coeff) return ta->coeff < tb->coeff ? -1 : 1;
    return 0;
}

static void die(const char *msg, const char *ctx) {
    fprintf(stderr, "PolyNormalize: %s%s%s\n",
            msg, ctx ? " near '" : "", ctx ? ctx : "");
    exit(2);
}

static int skip_ws(const char **p) {
    while (**p && isspace((unsigned char)**p)) (*p)++;
    return **p != 0;
}

static int parse_int(const char **p, long long *out) {
    const char *s = *p;
    int sign = 1;
    if (*s == '-') { sign = -1; s++; }
    else if (*s == '+') { s++; }
    if (!isdigit((unsigned char)*s)) return 0;
    long long v = 0;
    while (isdigit((unsigned char)*s)) {
        v = 10 * v + (*s - '0');
        s++;
    }
    *out = sign * v;
    *p = s;
    return 1;
}

/* Parse one factor: either "q" [ "^" NUMBER | "^{" NUMBER "}" ]
 *                  or "t" [ "^" NUMBER | "^{" NUMBER "}" ]
 * Returns 1 for q, 0 for t, -1 on error. exp set to 1 if no exponent. */
static int parse_factor(const char **p, long long *exp) {
    const char *s = *p;
    int kind;
    if      (*s == 'q') kind = 1;
    else if (*s == 't') kind = 0;
    else return -1;
    s++;
    long long e = 1;
    if (*s == '^') {
        s++;
        if (*s == '{') {
            s++;
            if (!parse_int(&s, &e)) return -1;
            if (*s != '}') return -1;
            s++;
        } else {
            if (!parse_int(&s, &e)) return -1;
        }
    }
    *exp = e;
    *p = s;
    return kind;
}

/* Parse one monomial term (optionally signed, with leading +/-).
 * Accepts: optional " + " or " - "; optional integer coefficient with
 * optional trailing '*'; then "q…" and/or "t…" in either order (missing
 * variable implies exponent 0 on the *output* side).
 * Sets have_q / have_t so a bare coefficient-only term like "7" can be
 * reported as has_q=0 has_t=0 and normalized to 7*q^0*t^0. */
static int parse_term(const char **p, Term *t, int *sign_was_set) {
    const char *s = *p;
    skip_ws(&s);
    long long sign = 1;
    *sign_was_set = 0;
    if (*s == '+' || *s == '-') {
        sign = (*s == '-') ? -1 : 1;
        *sign_was_set = 1;
        s++;
        skip_ws(&s);
    }
    long long coeff = sign;
    int saw_coeff = 0;
    /* Peek integer coefficient:
     *   "2q^3t" -> coeff 2  ; "2*q^3" -> coeff 2 ; "q^3" -> coeff sign*1 */
    const char *s2 = s;
    long long ctmp;
    if (isdigit((unsigned char)*s2)) {
        if (!parse_int(&s2, &ctmp)) return 0;
        coeff = sign * ctmp;
        saw_coeff = 1;
        if (*s2 == '*') s2++;
        s = s2;
        skip_ws(&s);
    }
    long long qexp = 0, texp = 0;
    int have_q = 0, have_t = 0;
    for (int k = 0; k < 2; k++) {
        long long e;
        int kd = parse_factor(&s, &e);
        if (kd == 1) { qexp = e; have_q = 1; }
        else if (kd == 0) { texp = e; have_t = 1; }
        else break;
        if (*s == '*') s++;
        skip_ws(&s);
    }
    if (!have_q && !have_t && !saw_coeff) return 0;
    t->coeff = coeff;
    t->q = qexp;
    t->t = texp;
    *p = s;
    return 1;
}

static void print_term(FILE *out, const Term *t, int first) {
    if (!first) {
        if (t->coeff < 0) fputc('-', out);
        else              fputc('+', out);
    } else {
        if (t->coeff < 0) fputc('-', out);
    }
    long long ac = t->coeff < 0 ? -t->coeff : t->coeff;
    if (ac != 1) fprintf(out, "%lld*", ac);
    fprintf(out, "q^%lld*t^%lld", t->q, t->t);
}

/* Combine duplicate terms (same q,t exponents) summing their coefficients. */
static size_t dedupe(Term *terms, size_t n) {
    if (n == 0) return 0;
    qsort(terms, n, sizeof(Term), cmp_terms);
    size_t w = 0;
    for (size_t i = 0; i < n; ) {
        size_t j = i;
        long long sum = 0;
        while (j < n && terms[j].q == terms[i].q && terms[j].t == terms[i].t) {
            sum += terms[j].coeff;
            j++;
        }
        if (sum != 0) {
            terms[w] = terms[i];
            terms[w].coeff = sum;
            w++;
        }
        i = j;
    }
    return w;
}

int main(int argc, char **argv) {
    const char *in_path = NULL, *out_path = NULL;
    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            printf("Usage: %s [input.poly] [output.poly]\n"
                   "       stdin/stdout if paths omitted.\n"
                   "Converts scanner-poly format -> gold 10a CSV format.\n", argv[0]);
            return 0;
        }
        else if (!in_path) in_path = argv[i];
        else if (!out_path) out_path = argv[i];
        else {
            fprintf(stderr, "extra argument: %s\n", argv[i]); return 2;
        }
    }
    FILE *in  = stdin;
    FILE *out = stdout;
    if (in_path  && !(in  = fopen(in_path,  "r"))) { perror(in_path);  return 2; }
    if (out_path && !(out = fopen(out_path, "w"))) { perror(out_path); return 2; }

    char *line = NULL;
    size_t cap = 0;
    ssize_t n;
    while ((n = getline(&line, &cap, in)) != -1) {
        /* Trim trailing newline(s) */
        while (n > 0 && (line[n-1] == '\n' || line[n-1] == '\r')) line[--n] = 0;
        if (n == 0) { fputc('\n', out); continue; }

        const char *p = line;
        skip_ws(&p);
        if (*p == 0) { fputc('\n', out); continue; }

        Term terms[2048];
        size_t nterms = 0;

        int sign_was_set = 0;
        Term tmp;
        /* First term: optional leading sign. Subsequent terms require a sign. */
        if (!parse_term(&p, &tmp, &sign_was_set)) {
            die("cannot parse first term of polynomial", line);
        }
        if (tmp.coeff != 0) terms[nterms++] = tmp;
        first = 0;
        skip_ws(&p);
        while (*p) {
            /* Must see + or - next. */
            if (*p != '+' && *p != '-') {
                die("expected + or - separator in polynomial", p);
            }
            if (!parse_term(&p, &tmp, &sign_was_set)) {
                die("cannot parse term", p);
            }
            if (tmp.coeff != 0) terms[nterms++] = tmp;
            skip_ws(&p);
        }
        if (nterms == 0) {
            fprintf(out, "0\n");
            continue;
        }
        size_t uniq = dedupe(terms, nterms);
        if (uniq == 0) { fprintf(out, "0\n"); continue; }
        for (size_t i = 0; i < uniq; i++) {
            print_term(out, &terms[i], (int)(i == 0));
        }
        fputc('\n', out);
    }
    free(line);
    if (in  != stdin)  fclose(in);
    if (out != stdout) fclose(out);
    return 0;
}
