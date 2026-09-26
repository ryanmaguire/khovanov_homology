#include "DTtoPD.h"
#include "PDScanner.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct DTPlanaritySearch {
  int n;
  int encounter_count; /* 2n */
  int dart_count;      /* 4n */

  int *lower_occurrence;
  int *upper_occurrence;
  int *signs;

  int *alpha;          /* edge involution on darts */
  int *sigma;          /* cyclic rotation at crossings */
  unsigned char *seen; /* scratch space for face counting */
} DTPlanaritySearch;

static void set_reason(char *reason, size_t reason_size, const char *text) {
  if (reason != NULL && reason_size > 0) {
    snprintf(reason, reason_size, "%s", text != NULL ? text : "DT conversion failed.");
  }
}

static int next_edge(int edge, int edge_count) {
  return (edge + 1 == edge_count) ? 0 : edge + 1;
}

static bool validate_numeric_dt(const int *dt, int n,
                                char *reason, size_t reason_size) {
  if (n < 0) {
    set_reason(reason, reason_size, "DT crossing count cannot be negative.");
    return false;
  }
  if (n == 0)
    return true;
  if (dt == NULL) {
    set_reason(reason, reason_size, "DT code is NULL.");
    return false;
  }
  if (n > DTTOPD_MAX_CROSSINGS) {
    if (reason != NULL && reason_size > 0) {
      snprintf(reason, reason_size,
               "DT-to-PD reconstruction currently supports at most %d crossings; got %d.",
               DTTOPD_MAX_CROSSINGS, n);
    }
    return false;
  }

  bool seen[DTTOPD_MAX_CROSSINGS + 1];
  memset(seen, 0, sizeof(seen));

  for (int i = 0; i < n; i++) {
    int value = dt[i];
    if (value == 0 || (value % 2) != 0) {
      if (reason != NULL && reason_size > 0) {
        snprintf(reason, reason_size,
                 "Invalid DT entry %d at position %d: entries must be nonzero even integers.",
                 value, i);
      }
      return false;
    }

    if (value < -2 * n || value > 2 * n) {
      if (reason != NULL && reason_size > 0) {
        snprintf(reason, reason_size,
                 "Invalid DT entry %d at position %d: absolute value must lie in 2..%d.",
                 value, i, 2 * n);
      }
      return false;
    }

    int abs_value = value < 0 ? -value : value;
    if (abs_value < 2) {
      if (reason != NULL && reason_size > 0) {
        snprintf(reason, reason_size,
                 "Invalid DT entry %d at position %d: absolute value must lie in 2..%d.",
                 value, i, 2 * n);
      }
      return false;
    }

    int even_index = abs_value / 2;
    if (seen[even_index]) {
      if (reason != NULL && reason_size > 0) {
        snprintf(reason, reason_size,
                 "Invalid DT code: absolute even label %d occurs more than once.",
                 abs_value);
      }
      return false;
    }
    seen[even_index] = true;
  }

  for (int k = 1; k <= n; k++) {
    if (!seen[k]) {
      if (reason != NULL && reason_size > 0) {
        snprintf(reason, reason_size,
                 "Invalid DT code: even label %d is missing.", 2 * k);
      }
      return false;
    }
  }

  return true;
}

bool DTtoPD_decodeAlphabetical(const char *code,
                               int **dt_out,
                               int *crossing_count_out,
                               char *reason,
                               size_t reason_size) {
  if (dt_out == NULL || crossing_count_out == NULL) {
    set_reason(reason, reason_size, "DT decoder requires non-NULL output pointers.");
    return false;
  }

  *dt_out = NULL;
  *crossing_count_out = 0;

  if (code == NULL) {
    set_reason(reason, reason_size, "Alphabetical DT code is NULL.");
    return false;
  }

  int count = 0;
  for (const unsigned char *p = (const unsigned char *)code; *p != '\0'; p++) {
    if (isspace(*p))
      continue;
    if (!((*p >= 'a' && *p <= 'z') || (*p >= 'A' && *p <= 'Z'))) {
      if (reason != NULL && reason_size > 0) {
        snprintf(reason, reason_size,
                 "Invalid alphabetical DT character '%c'; expected ASCII letters only.", *p);
      }
      return false;
    }
    count++;
  }

  if (count > DTTOPD_MAX_CROSSINGS) {
    if (reason != NULL && reason_size > 0) {
      snprintf(reason, reason_size,
               "Alphabetical DT code has %d crossings; current DT-to-PD limit is %d.",
               count, DTTOPD_MAX_CROSSINGS);
    }
    return false;
  }

  if (count == 0)
    return true;

  int *dt = (int *)malloc((size_t)count * sizeof(int));
  if (dt == NULL) {
    set_reason(reason, reason_size, "Out of memory while decoding alphabetical DT code.");
    return false;
  }

  int out = 0;
  for (const unsigned char *p = (const unsigned char *)code; *p != '\0'; p++) {
    if (isspace(*p))
      continue;

    if (*p >= 'a' && *p <= 'z') {
      dt[out++] = 2 * ((int)(*p - 'a') + 1);
    } else if (*p >= 'A' && *p <= 'Z') {
      dt[out++] = -2 * ((int)(*p - 'A') + 1);
    } else {
      free(dt);
      set_reason(reason, reason_size,
                 "Alphabetical DT input must use ASCII letters A-Z or a-z.");
      return false;
    }
  }

  if (!validate_numeric_dt(dt, count, reason, reason_size)) {
    free(dt);
    return false;
  }

  *dt_out = dt;
  *crossing_count_out = count;
  return true;
}

/*
 * Darts are encoded as 2*occurrence + side, where side 0 is the incoming
 * half-edge and side 1 is the outgoing half-edge along the oriented knot.
 *
 * alpha pairs the two darts belonging to the same diagram edge.
 */
static void build_alpha(DTPlanaritySearch *s) {
  int m = s->encounter_count;
  for (int occurrence = 0; occurrence < m; occurrence++) {
    int prev = occurrence == 0 ? m - 1 : occurrence - 1;
    int next = occurrence + 1 == m ? 0 : occurrence + 1;

    s->alpha[2 * occurrence] = 2 * prev + 1;      /* incoming <-> previous outgoing */
    s->alpha[2 * occurrence + 1] = 2 * next;      /* outgoing <-> next incoming */
  }
}

static void set_rotation_cycle(int *sigma, int a, int b, int c, int d) {
  sigma[a] = b;
  sigma[b] = c;
  sigma[c] = d;
  sigma[d] = a;
}

/*
 * Build the orientable rotation system for the current crossing signs.
 *
 * For sign +1 the cyclic order, starting at the incoming lower strand, is
 *     lower_in, upper_out, lower_out, upper_in.
 * For sign -1 it is
 *     lower_in, upper_in, lower_out, upper_out.
 *
 * These are exactly the two local rotations compatible with the two strands
 * remaining opposite at a crossing.
 */
static void build_sigma(DTPlanaritySearch *s) {
  for (int i = 0; i < s->n; i++) {
    int lower = s->lower_occurrence[i];
    int upper = s->upper_occurrence[i];

    int lower_in = 2 * lower;
    int lower_out = lower_in + 1;
    int upper_in = 2 * upper;
    int upper_out = upper_in + 1;

    if (s->signs[i] > 0) {
      set_rotation_cycle(s->sigma, lower_in, upper_out, lower_out, upper_in);
    } else {
      set_rotation_cycle(s->sigma, lower_in, upper_in, lower_out, upper_out);
    }
  }
}

/* Count faces of the orientable ribbon graph defined by sigma and alpha. */
static int count_faces(DTPlanaritySearch *s) {
  build_sigma(s);
  memset(s->seen, 0, (size_t)s->dart_count * sizeof(unsigned char));

  int faces = 0;
  for (int dart = 0; dart < s->dart_count; dart++) {
    if (s->seen[dart])
      continue;

    faces++;
    int current = dart;
    while (!s->seen[current]) {
      s->seen[current] = 1;
      current = s->sigma[s->alpha[current]];
    }
  }
  return faces;
}

/*
 * A connected 4-valent n-crossing knot shadow has V=n and E=2n.
 * A rotation system is planar exactly when Euler's formula gives
 *
 *     V - E + F = 2  <=>  F = n + 2.
 *
 * DT notation does not distinguish the global mirror.  Negating every local
 * crossing sign reverses every crossing rotation and preserves planarity, so
 * we fix crossing 0 to +1 and search only the remaining n-1 choices.
 */
static bool search_planar_rotation(DTPlanaritySearch *s, int crossing) {
  if (crossing == s->n)
    return count_faces(s) == s->n + 2;

  s->signs[crossing] = 1;
  if (search_planar_rotation(s, crossing + 1))
    return true;

  s->signs[crossing] = -1;
  if (search_planar_rotation(s, crossing + 1))
    return true;

  return false;
}

static bool find_planar_signs(const int *dt, int n, int *signs,
                              char *reason, size_t reason_size) {
  if (n == 0)
    return true;

  DTPlanaritySearch s;
  memset(&s, 0, sizeof(s));
  s.n = n;
  s.encounter_count = 2 * n;
  s.dart_count = 4 * n;
  s.signs = signs;

  s.lower_occurrence = (int *)malloc((size_t)n * sizeof(int));
  s.upper_occurrence = (int *)malloc((size_t)n * sizeof(int));
  s.alpha = (int *)malloc((size_t)s.dart_count * sizeof(int));
  s.sigma = (int *)malloc((size_t)s.dart_count * sizeof(int));
  s.seen = (unsigned char *)malloc((size_t)s.dart_count * sizeof(unsigned char));

  if (s.lower_occurrence == NULL || s.upper_occurrence == NULL ||
      s.alpha == NULL || s.sigma == NULL || s.seen == NULL) {
    free(s.lower_occurrence);
    free(s.upper_occurrence);
    free(s.alpha);
    free(s.sigma);
    free(s.seen);
    set_reason(reason, reason_size, "Out of memory while reconstructing DT planar embedding.");
    return false;
  }

  for (int i = 0; i < n; i++) {
    int odd_occurrence = 2 * i;            /* DT odd label 2i+1, converted to 0-based */
    int even_occurrence = abs(dt[i]) - 1;  /* signed even DT label, converted to 0-based */

    /*
     * DT sign convention: a negative even label means that the even-numbered
     * encounter is the overpass.  Therefore the lower (under) encounter is the
     * even one for dt>0, and the odd one for dt<0.
     */
    if (dt[i] > 0) {
      s.lower_occurrence[i] = even_occurrence;
      s.upper_occurrence[i] = odd_occurrence;
    } else {
      s.lower_occurrence[i] = odd_occurrence;
      s.upper_occurrence[i] = even_occurrence;
    }
  }

  build_alpha(&s);

  signs[0] = 1; /* deterministic choice between the two global mirrors */
  bool ok = search_planar_rotation(&s, 1);

  free(s.lower_occurrence);
  free(s.upper_occurrence);
  free(s.alpha);
  free(s.sigma);
  free(s.seen);

  if (!ok) {
    set_reason(reason, reason_size,
               "DT code did not admit a planar rotation system under the supported classical-knot convention.");
  }
  return ok;
}

static void fill_pd_crossing(int out[4], int lower, int upper,
                             int crossing_sign, int edge_count) {
  int lower_in = lower;
  int lower_out = next_edge(lower, edge_count);
  int upper_in = upper;
  int upper_out = next_edge(upper, edge_count);

  out[0] = lower_in;
  out[2] = lower_out;

  if (crossing_sign > 0) {
    out[1] = upper_out;
    out[3] = upper_in;
  } else {
    out[1] = upper_in;
    out[3] = upper_out;
  }
}

static bool validate_built_pd(const PDDiagram *diagram,
                              char *reason, size_t reason_size) {
  if (diagram->crossing_count == 0)
    return true;

  int edge_count = diagram->edge_count;
  int *uses = (int *)calloc((size_t)edge_count, sizeof(int));
  if (uses == NULL) {
    set_reason(reason, reason_size, "Out of memory while validating reconstructed PD.");
    return false;
  }

  for (int i = 0; i < diagram->crossing_count; i++) {
    if (diagram->signs[i] != 1 && diagram->signs[i] != -1) {
      free(uses);
      set_reason(reason, reason_size, "Internal DT-to-PD error: invalid reconstructed crossing sign.");
      return false;
    }

    for (int j = 0; j < 4; j++) {
      int edge = diagram->crossings[i][j];
      if (edge < 0 || edge >= edge_count) {
        free(uses);
        set_reason(reason, reason_size, "Internal DT-to-PD error: reconstructed edge label is out of range.");
        return false;
      }
      uses[edge]++;
    }
  }

  for (int edge = 0; edge < edge_count; edge++) {
    if (uses[edge] != 2) {
      int count = uses[edge];
      free(uses);
      if (reason != NULL && reason_size > 0) {
        snprintf(reason, reason_size,
                 "Internal DT-to-PD error: edge %d occurs %d times instead of twice.",
                 edge, count);
      }
      return false;
    }
  }

  free(uses);
  return true;
}

bool DTtoPD_fromNumeric(PDDiagram *diagram,
                        const int *dt,
                        int crossing_count,
                        char *reason,
                        size_t reason_size) {
  if (diagram == NULL) {
    set_reason(reason, reason_size, "DT-to-PD output diagram is NULL.");
    return false;
  }

  memset(diagram, 0, sizeof(*diagram));

  if (!validate_numeric_dt(dt, crossing_count, reason, reason_size))
    return false;

  if (crossing_count == 0)
    return true;

  diagram->crossing_count = crossing_count;
  diagram->edge_count = 2 * crossing_count;
  diagram->crossings =
      (int (*)[4])malloc((size_t)crossing_count * sizeof(*diagram->crossings));
  diagram->signs = (int *)malloc((size_t)crossing_count * sizeof(int));

  if (diagram->crossings == NULL || diagram->signs == NULL) {
    PDDiagram_free(diagram);
    set_reason(reason, reason_size, "Out of memory while allocating reconstructed PD.");
    return false;
  }

  if (!find_planar_signs(dt, crossing_count, diagram->signs,
                         reason, reason_size)) {
    PDDiagram_free(diagram);
    return false;
  }

  for (int i = 0; i < crossing_count; i++) {
    int odd_occurrence = 2 * i;
    int even_occurrence = abs(dt[i]) - 1;
    int lower = dt[i] > 0 ? even_occurrence : odd_occurrence;
    int upper = dt[i] > 0 ? odd_occurrence : even_occurrence;

    fill_pd_crossing(diagram->crossings[i], lower, upper,
                     diagram->signs[i], diagram->edge_count);
  }

  if (!validate_built_pd(diagram, reason, reason_size)) {
    PDDiagram_free(diagram);
    return false;
  }

  return true;
}

bool DTtoPD_fromAlphabetical(PDDiagram *diagram,
                             const char *code,
                             char *reason,
                             size_t reason_size) {
  int *dt = NULL;
  int crossings = 0;

  if (!DTtoPD_decodeAlphabetical(code, &dt, &crossings,
                                 reason, reason_size)) {
    return false;
  }

  bool ok = DTtoPD_fromNumeric(diagram, dt, crossings, reason, reason_size);
  free(dt);
  return ok;
}
