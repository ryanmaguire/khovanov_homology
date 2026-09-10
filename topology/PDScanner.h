#ifndef PDSCANNER_H
#define PDSCANNER_H

#include "Komplex.h"
#include <stdbool.h>
#include <stddef.h>

typedef struct PDDiagram {
  int crossing_count;
  int edge_count;
  int (*crossings)[4];
  int *signs;
} PDDiagram;

bool PDDiagram_parse(PDDiagram *diagram, const char *pd_text,
                     const char *sign_text, bool infer_javakh_signs,
                     char *reason, size_t reason_size);
void PDDiagram_free(PDDiagram *diagram);

Komplex *PDScanner_build(const PDDiagram *diagram, bool reorder_crossings,
                         bool verify_d_squared, char *reason,
                         size_t reason_size);

#endif
