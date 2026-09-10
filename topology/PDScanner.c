#include "PDScanner.h"
#include "TangleCompose.h"
#include "CannedCobordismImpl.h"
#include "LCCC.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void reduce_complex(Komplex *k) {
  if (k == NULL) return;
  Komplex_deloop(k);
  for (int i = 0; i < k->length - 1; i++)
    if (k->differentials[i] != NULL) CobMatrix_reduce(k->differentials[i]);
  Komplex_greedyReduce(k);
  for (int i = 0; i < k->length - 1; i++)
    if (k->differentials[i] != NULL) CobMatrix_reduce(k->differentials[i]);
}

static bool parse_integer_list(const char *text, int **values, int *count) {
  *values = NULL;
  *count = 0;
  if (text == NULL) return true;

  int capacity = 16;
  int *out = (int *)malloc((size_t)capacity * sizeof(int));
  if (out == NULL) return false;

  const char *p = text;
  while (*p != '\0') {
    while (*p != '\0' && !isdigit((unsigned char)*p) && *p != '-' && *p != '+') p++;
    if (*p == '\0') break;

    char *end = NULL;
    long value = strtol(p, &end, 10);
    if (end == p) {
      p++;
      continue;
    }
    if (value < -2147483647L - 1L || value > 2147483647L) {
      free(out);
      return false;
    }
    if (*count == capacity) {
      capacity *= 2;
      int *grown = (int *)realloc(out, (size_t)capacity * sizeof(int));
      if (grown == NULL) {
        free(out);
        return false;
      }
      out = grown;
    }
    out[(*count)++] = (int)value;
    p = end;
  }

  *values = out;
  return true;
}

static bool infer_signs(const int (*raw)[4], int crossings, int *signs,
                        char *reason, size_t reason_size) {
  for (int i = 0; i < crossings; i++) {
    int d13 = raw[i][1] - raw[i][3];
    int d31 = raw[i][3] - raw[i][1];
    if (d13 == 1 || d31 > 1)
      signs[i] = 1;
    else if (d31 == 1 || d13 > 1)
      signs[i] = -1;
    else {
      snprintf(reason, reason_size,
               "JavaKh sign inference failed at crossing %d; provide --signs explicitly.",
               i);
      return false;
    }
  }
  return true;
}

bool PDDiagram_parse(PDDiagram *diagram, const char *pd_text,
                     const char *sign_text, bool infer_javakh_signs,
                     char *reason, size_t reason_size) {
  if (diagram == NULL) return false;
  memset(diagram, 0, sizeof(*diagram));

  int *flat = NULL;
  int flat_count = 0;
  if (!parse_integer_list(pd_text, &flat, &flat_count)) {
    snprintf(reason, reason_size, "Failed to parse PD integers.");
    return false;
  }
  if ((flat_count % 4) != 0) {
    free(flat);
    snprintf(reason, reason_size,
             "A PD must contain four edge labels per crossing; parsed %d integers.",
             flat_count);
    return false;
  }

  diagram->crossing_count = flat_count / 4;
  if (diagram->crossing_count == 0) {
    free(flat);
    return true;
  }

  diagram->crossings =
      (int (*)[4])malloc((size_t)diagram->crossing_count * sizeof(*diagram->crossings));
  diagram->signs = (int *)malloc((size_t)diagram->crossing_count * sizeof(int));
  if (diagram->crossings == NULL || diagram->signs == NULL) {
    free(flat);
    PDDiagram_free(diagram);
    snprintf(reason, reason_size, "Out of memory while parsing the PD.");
    return false;
  }

  for (int i = 0; i < diagram->crossing_count; i++)
    for (int j = 0; j < 4; j++) diagram->crossings[i][j] = flat[4 * i + j];
  free(flat);

  if (sign_text != NULL) {
    int *parsed_signs = NULL;
    int sign_count = 0;
    if (!parse_integer_list(sign_text, &parsed_signs, &sign_count) ||
        sign_count != diagram->crossing_count) {
      free(parsed_signs);
      PDDiagram_free(diagram);
      snprintf(reason, reason_size,
               "--signs must contain exactly one +1 or -1 per crossing.");
      return false;
    }
    for (int i = 0; i < sign_count; i++) {
      if (parsed_signs[i] != 1 && parsed_signs[i] != -1) {
        free(parsed_signs);
        PDDiagram_free(diagram);
        snprintf(reason, reason_size, "Crossing signs must be +1 or -1.");
        return false;
      }
      diagram->signs[i] = parsed_signs[i];
    }
    free(parsed_signs);
  } else if (infer_javakh_signs) {
    if (!infer_signs((const int (*)[4])diagram->crossings,
                     diagram->crossing_count, diagram->signs, reason,
                     reason_size)) {
      PDDiagram_free(diagram);
      return false;
    }
  } else {
    PDDiagram_free(diagram);
    snprintf(reason, reason_size,
             "PD mode requires --signs unless --javakh-signs is specified.");
    return false;
  }

  int max_unique = 2 * diagram->crossing_count;
  int *labels = (int *)malloc((size_t)max_unique * sizeof(int));
  int *counts = (int *)calloc((size_t)max_unique, sizeof(int));
  if (labels == NULL || counts == NULL) {
    free(labels);
    free(counts);
    PDDiagram_free(diagram);
    snprintf(reason, reason_size, "Out of memory while normalizing PD edges.");
    return false;
  }

  int unique = 0;
  for (int i = 0; i < diagram->crossing_count; i++) {
    for (int j = 0; j < 4; j++) {
      int label = diagram->crossings[i][j];
      int dense = -1;
      for (int k = 0; k < unique; k++)
        if (labels[k] == label) { dense = k; break; }
      if (dense < 0) {
        if (unique >= max_unique) {
          free(labels);
          free(counts);
          PDDiagram_free(diagram);
          snprintf(reason, reason_size,
                   "A closed PD has more distinct edge labels than expected.");
          return false;
        }
        dense = unique;
        labels[unique++] = label;
      }
      counts[dense]++;
      diagram->crossings[i][j] = dense;
    }
  }

  for (int i = 0; i < unique; i++) {
    if (counts[i] != 2) {
      int label = labels[i];
      int count = counts[i];
      free(labels);
      free(counts);
      PDDiagram_free(diagram);
      snprintf(reason, reason_size,
               "PD edge label %d occurs %d times; every closed-PD edge must occur exactly twice.",
               label, count);
      return false;
    }
  }
  diagram->edge_count = unique;
  free(labels);
  free(counts);
  return true;
}

void PDDiagram_free(PDDiagram *diagram) {
  if (diagram == NULL) return;
  free(diagram->crossings);
  free(diagram->signs);
  memset(diagram, 0, sizeof(*diagram));
}

static Komplex *make_local_crossing(int sign) {
  Komplex *k = Komplex_create(2);
  if (k == NULL) return NULL;

  Cap *zero = Cap_create(4, 0);
  Cap *one = Cap_create(4, 0);
  if (zero == NULL || one == NULL) return NULL;
  zero->pairings[0] = 1; zero->pairings[1] = 0;
  zero->pairings[2] = 3; zero->pairings[3] = 2;
  one->pairings[0] = 3; one->pairings[3] = 0;
  one->pairings[1] = 2; one->pairings[2] = 1;

  SmoothingColumn *c0 = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  SmoothingColumn *c1 = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  if (c0 == NULL || c1 == NULL) return NULL;
  c0->n = c1->n = 1;
  c0->numbers = (int *)malloc(sizeof(int));
  c1->numbers = (int *)malloc(sizeof(int));
  c0->smoothings = (Cap **)malloc(sizeof(Cap *));
  c1->smoothings = (Cap **)malloc(sizeof(Cap *));
  if (c0->numbers == NULL || c1->numbers == NULL ||
      c0->smoothings == NULL || c1->smoothings == NULL) return NULL;

  c0->smoothings[0] = zero;
  c1->smoothings[0] = one;
  c0->numbers[0] = sign > 0 ? 0 : -1;
  c1->numbers[0] = sign > 0 ? 1 : 0;
  k->chain_groups[0] = c0;
  k->chain_groups[1] = c1;

  CannedCobordismImplData *impl = CannedCobordismImpl_create(zero, one);
  if (impl == NULL) return NULL;
  impl->ncc = impl->nbc;
  for (int i = 0; i < impl->nbc; i++) impl->connectedComponent[i] = i;
  impl->dots = (int *)calloc((size_t)impl->ncc, sizeof(int));
  impl->genus = (int *)calloc((size_t)impl->ncc, sizeof(int));
  if (impl->ncc > 0 && (impl->dots == NULL || impl->genus == NULL)) return NULL;

  CannedCobordism *saddle = CannedCobordismImpl_as_CannedCobordism(impl);
  CobMatrix *d = CobMatrix_create(c0, c1, true);
  if (saddle == NULL || d == NULL) return NULL;
  CobMatrix_putEntry(d, 0, 0, LCCC_createSingle(saddle, 1));
  k->differentials[0] = d;
  return k;
}

static Komplex *make_arc_complex(void) {
  Komplex *k = Komplex_create(1);
  if (k == NULL) return NULL;
  SmoothingColumn *c = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  if (c == NULL) return NULL;
  c->n = 1;
  c->numbers = (int *)malloc(sizeof(int));
  c->smoothings = (Cap **)malloc(sizeof(Cap *));
  if (c->numbers == NULL || c->smoothings == NULL) return NULL;
  c->numbers[0] = 0;
  c->smoothings[0] = Cap_create(2, 0);
  if (c->smoothings[0] == NULL) return NULL;
  c->smoothings[0]->pairings[0] = 1;
  c->smoothings[0]->pairings[1] = 0;
  k->chain_groups[0] = c;
  return k;
}

static Komplex *make_unknot_complex(void) {
  Komplex *k = Komplex_create(1);
  if (k == NULL) return NULL;
  SmoothingColumn *c = (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  if (c == NULL) return NULL;
  c->n = 1;
  c->numbers = (int *)malloc(sizeof(int));
  c->smoothings = (Cap **)malloc(sizeof(Cap *));
  if (c->numbers == NULL || c->smoothings == NULL) return NULL;
  c->numbers[0] = 0;
  c->smoothings[0] = Cap_create(0, 1);
  if (c->smoothings[0] == NULL) return NULL;
  k->chain_groups[0] = c;
  return k;
}

static bool attachment(const int *frontier, int frontier_count,
                       const int crossing[4], const bool *seen,
                       int *join_count, int *frontier_start,
                       int *crossing_start) {
  int ncon = 0;
  for (int j = 0; j < 4; j++) if (seen[crossing[j]]) ncon++;
  if (frontier_count == 0) {
    if (ncon != 0) return false;
    *join_count = 0;
    *frontier_start = 0;
    *crossing_start = 0;
    return true;
  }
  if (ncon == 0 || ncon > frontier_count) return false;

  for (int start = 0; start < frontier_count; start++) {
    for (int ks = 0; ks < 4; ks++) {
      bool good = true;
      for (int j = 0; j < ncon; j++) {
        if (crossing[(ks - j + 8) % 4] !=
            frontier[(start + j) % frontier_count]) {
          good = false;
          break;
        }
      }
      if (good) {
        *join_count = ncon;
        *frontier_start = start;
        *crossing_start = (ks + 4 - ncon + 1) % 4;
        return true;
      }
    }
  }
  return false;
}

static int choose_crossing(const PDDiagram *diagram, const int *frontier,
                           int frontier_count, const bool *seen,
                           const bool *done, bool reorder, int *join_count,
                           int *frontier_start, int *crossing_start) {
  int best = -1;
  int best_join = -1;
  for (int i = 0; i < diagram->crossing_count; i++) {
    if (done[i]) continue;
    int jc = 0, fs = 0, cs = 0;
    if (!attachment(frontier, frontier_count, diagram->crossings[i], seen,
                    &jc, &fs, &cs)) continue;
    if (!reorder) {
      *join_count = jc;
      *frontier_start = fs;
      *crossing_start = cs;
      return i;
    }
    if (jc > best_join) {
      best = i;
      best_join = jc;
      *join_count = jc;
      *frontier_start = fs;
      *crossing_start = cs;
    }
  }
  return best;
}

static int *updated_frontier(const int *frontier, int frontier_count,
                             int frontier_start, const int crossing[4],
                             int crossing_start, int join_count,
                             int *new_count) {
  *new_count = frontier_count + 4 - 2 * join_count;
  int *out = *new_count > 0 ? (int *)malloc((size_t)*new_count * sizeof(int)) : NULL;
  if (*new_count > 0 && out == NULL) return NULL;

  int n = 0;
  for (; n < frontier_count - join_count; n++)
    out[n] = frontier[(frontier_start + join_count + n) % frontier_count];
  for (int j = 0; j < 4 - join_count; j++, n++)
    out[n] = crossing[(crossing_start + join_count + j) % 4];
  return out;
}

static int find_adjacent_duplicate(const int *frontier, int count) {
  if (count < 2) return -1;
  for (int i = 0; i < count; i++)
    if (frontier[i] == frontier[(i + 1) % count]) return i;
  return -1;
}

static bool has_duplicate(const int *frontier, int count) {
  for (int i = 0; i < count; i++)
    for (int j = i + 1; j < count; j++)
      if (frontier[i] == frontier[j]) return true;
  return false;
}

static bool remove_cyclic_pair(int *frontier, int *count, int start) {
  int old = *count;
  int *tmp = old > 2 ? (int *)malloc((size_t)(old - 2) * sizeof(int)) : NULL;
  if (old > 2 && tmp == NULL) return false;
  int n = 0;
  for (int k = 2; k < old; k++) tmp[n++] = frontier[(start + k) % old];
  for (int i = 0; i < n; i++) frontier[i] = tmp[i];
  free(tmp);
  *count = old - 2;
  return true;
}

static bool check_boundary(const Komplex *k, int expected, char *reason,
                           size_t reason_size) {
  for (int h = 0; h < k->length; h++) {
    SmoothingColumn *col = k->chain_groups[h];
    if (col == NULL) continue;
    for (int i = 0; i < col->n; i++) {
      if (col->smoothings[i] == NULL || col->smoothings[i]->n != expected) {
        snprintf(reason, reason_size,
                 "PD frontier has %d edges but C_%d[%d] has boundary size %d.",
                 expected, h, i,
                 col->smoothings[i] == NULL ? -1 : col->smoothings[i]->n);
        return false;
      }
    }
  }
  return true;
}

static bool close_internal_edges(Komplex **complex, int *frontier,
                                 int *frontier_count, bool verify_d_squared,
                                 char *reason, size_t reason_size) {
  while (has_duplicate(frontier, *frontier_count)) {
    int start = find_adjacent_duplicate(frontier, *frontier_count);
    if (start < 0) {
      snprintf(reason, reason_size,
               "The processed PD frontier contains a non-adjacent internal edge; the current region is not represented by one simply connected frontier.");
      return false;
    }

    Komplex *arc = make_arc_complex();
    Komplex *next =
        Komplex_compose_partial_tangles(*complex, start, arc, 0, 2);
    Komplex_free(*complex);
    Komplex_free(arc);
    if (next == NULL) {
      snprintf(reason, reason_size, "Failed to close an internal PD edge.");
      return false;
    }
    *complex = next;
    if (!remove_cyclic_pair(frontier, frontier_count, start)) {
      snprintf(reason, reason_size,
               "Out of memory while removing a closed PD frontier edge.");
      return false;
    }
    reduce_complex(*complex);
    if (verify_d_squared && !Komplex_verify_d_squared(*complex)) {
      snprintf(reason, reason_size,
               "d^2 != 0 after closing an internal PD edge.");
      return false;
    }
  }
  return true;
}

Komplex *PDScanner_build(const PDDiagram *diagram, bool reorder_crossings,
                         bool verify_d_squared, char *reason,
                         size_t reason_size) {
  if (diagram == NULL) return NULL;
  if (diagram->crossing_count == 0) return make_unknot_complex();

  bool *seen = (bool *)calloc((size_t)diagram->edge_count, sizeof(bool));
  bool *done = (bool *)calloc((size_t)diagram->crossing_count, sizeof(bool));
  if (seen == NULL || done == NULL) {
    free(seen); free(done);
    snprintf(reason, reason_size, "Out of memory while starting the PD scan.");
    return NULL;
  }

  int *frontier = NULL;
  int frontier_count = 0;
  Komplex *current = NULL;

  for (int step = 0; step < diagram->crossing_count; step++) {
    int join_count = 0, frontier_start = 0, crossing_start = 0;
    int best = choose_crossing(diagram, frontier, frontier_count, seen, done,
                               reorder_crossings, &join_count, &frontier_start,
                               &crossing_start);
    if (best < 0) {
      snprintf(reason, reason_size,
               "No unprocessed crossing attaches as a consecutive block to the current PD frontier.");
      free(frontier); free(seen); free(done); Komplex_free(current);
      return NULL;
    }

    Komplex *cross = make_local_crossing(diagram->signs[best]);
    if (cross == NULL) {
      snprintf(reason, reason_size, "Failed to construct crossing %d.", best);
      free(frontier); free(seen); free(done); Komplex_free(current);
      return NULL;
    }

    if (current == NULL) {
      current = cross;
    } else {
      Komplex *next = Komplex_compose_partial_tangles(
          current, frontier_start, cross, crossing_start, join_count);
      Komplex_free(current);
      Komplex_free(cross);
      if (next == NULL) {
        snprintf(reason, reason_size, "Partial composition failed at crossing %d.", best);
        free(frontier); free(seen); free(done);
        return NULL;
      }
      current = next;
    }

    int new_count = 0;
    int *new_frontier = updated_frontier(
        frontier, frontier_count, frontier_start, diagram->crossings[best],
        crossing_start, join_count, &new_count);
    if (new_count > 0 && new_frontier == NULL) {
      snprintf(reason, reason_size, "Out of memory while updating the PD frontier.");
      free(frontier); free(seen); free(done); Komplex_free(current);
      return NULL;
    }
    free(frontier);
    frontier = new_frontier;
    frontier_count = new_count;

    done[best] = true;
    for (int j = 0; j < 4; j++) seen[diagram->crossings[best][j]] = true;

    if (!close_internal_edges(&current, frontier, &frontier_count,
                              verify_d_squared, reason, reason_size)) {
      free(frontier); free(seen); free(done); Komplex_free(current);
      return NULL;
    }

    reduce_complex(current);
    if (!check_boundary(current, frontier_count, reason, reason_size)) {
      free(frontier); free(seen); free(done); Komplex_free(current);
      return NULL;
    }
    if (verify_d_squared && !Komplex_verify_d_squared(current)) {
      snprintf(reason, reason_size, "d^2 != 0 after scanning crossing %d.", best);
      free(frontier); free(seen); free(done); Komplex_free(current);
      return NULL;
    }
  }

  free(seen);
  free(done);
  if (frontier_count != 0) {
    snprintf(reason, reason_size,
             "PD scan ended with %d exposed frontier edges instead of a closed diagram.",
             frontier_count);
    free(frontier);
    Komplex_free(current);
    return NULL;
  }
  free(frontier);
  return current;
}