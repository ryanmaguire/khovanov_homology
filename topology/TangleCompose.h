#ifndef TANGLECOMPOSE_H
#define TANGLECOMPOSE_H

#include "Komplex.h"

Komplex *Komplex_compose_partial_tangles(Komplex *left, int left_start,
                                         Komplex *right, int right_start,
                                         int join_count);

#endif
