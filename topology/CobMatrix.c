#include "../KhMemory.h"
#include "CobMatrix.h"
#include "Cap.h"
#include <string.h>
#include <assert.h>
#include <stdlib.h>

SmoothingColumn *SmoothingColumn_create(void) {
  SmoothingColumn *col = kh_calloc(1, sizeof(*col));
  col->references = 1;
  return col;
}
SmoothingColumn *SmoothingColumn_retain(SmoothingColumn *col) {
  if (col) {
    if (col->references == SIZE_MAX) kh_fatal("column reference count overflow");
    ++col->references;
  }
  return col;
}
void SmoothingColumn_free(SmoothingColumn *col) {
  if (!col || --col->references) return;
  for (int i = 0; i < col->n; ++i) Cap_free(col->smoothings[i]);
  free(col->numbers);
  free(col->smoothings);
  free(col);
}

// Helper to create a new linked list node
static MatrixEntry *create_entry(int col_idx, LCCC *value) {
  MatrixEntry *entry = (MatrixEntry *)kh_malloc(sizeof(MatrixEntry));
  entry->column_index = col_idx;
  entry->value = value;
  entry->next = NULL;
  return entry;
}

// Helper to free a matrix row list.
//
// Matrix entries own their LCCC values. Row/column extraction moves entries
// rather than sharing them, and CobMatrix_add clones source values before
// inserting them into the destination.
static void free_matrix_row(MatrixRow *row) {
  MatrixEntry *current = row->head;
  while (current != NULL) {
    MatrixEntry *next = current->next;
    LCCC_free(current->value);
    free(current);
    current = next;
  }
  row->head = NULL;
}

// ----------------------------------------------------
// Constructor & Destructor
// ----------------------------------------------------

CobMatrix *CobMatrix_create(SmoothingColumn *source, SmoothingColumn *target,
                            bool shared) {
  CobMatrix *m = (CobMatrix *)kh_malloc(sizeof(CobMatrix));
  if (shared) {
    m->source = SmoothingColumn_retain(source);
    m->target = SmoothingColumn_retain(target);
  } else {
    m->source = SmoothingColumn_clone(source);
    m->target = SmoothingColumn_clone(target);
  }

  m->entries = (MatrixRow *)kh_calloc(m->target->n, sizeof(MatrixRow));
  for (int i = 0; i < m->target->n; i++) {
    m->entries[i].head = NULL;
  }

  return m;
}

void CobMatrix_free(CobMatrix *m) {
  if (!m)
    return;
  for (int i = 0; i < m->target->n; i++) {
    free_matrix_row(&m->entries[i]);
  }
  free(m->entries);
  SmoothingColumn_free(m->source);
  SmoothingColumn_free(m->target);
  free(m);
}

// ----------------------------------------------------
// Core Matrix Operations
// ----------------------------------------------------

void CobMatrix_putEntry(CobMatrix *m, int row_idx, int col_idx, LCCC *lc) {
  // Ownership of lc is transferred to the matrix.
  if (lc == NULL)
    return;
  if (LCCC_isZero(lc)) {
    LCCC_free(lc);
    return;
  }

  MatrixRow *row = &m->entries[row_idx];
  MatrixEntry *new_entry = create_entry(col_idx, lc);

  // Insert at the head. A later optimization may keep rows sorted, but sparse
  // linked rows remain correct without ordering.
  new_entry->next = row->head;
  row->head = new_entry;
}

void CobMatrix_addEntry(CobMatrix *m, int row_idx, int col_idx, LCCC *t) {
  // Ownership of t is transferred to this function whether the entry is
  // inserted, merged, or vanishes.
  if (t == NULL)
    return;
  if (LCCC_isZero(t)) {
    LCCC_free(t);
    return;
  }

  MatrixRow *row = &m->entries[row_idx];
  MatrixEntry *prev = NULL;
  MatrixEntry *current = row->head;

  while (current != NULL) {
    if (current->column_index == col_idx) {
      LCCC *sum = LCCC_add(current->value, t);
      LCCC_free(t);

      if (sum == NULL || LCCC_isZero(sum)) {
        if (prev) {
          prev->next = current->next;
        } else {
          row->head = current->next;
        }
        LCCC_free(current->value);
        LCCC_free(sum);
        free(current);
      } else {
        LCCC_free(current->value);
        current->value = sum;
      }
      return;
    }
    prev = current;
    current = current->next;
  }

  MatrixEntry *new_entry = create_entry(col_idx, t);
  new_entry->next = row->head;
  row->head = new_entry;
}

LCCC **CobMatrix_unpackRow(CobMatrix *m, int row_idx) {
  LCCC **unpacked = (LCCC **)kh_calloc(m->source->n, sizeof(LCCC *));
  MatrixEntry *current = m->entries[row_idx].head;
  while (current != NULL) {
    unpacked[current->column_index] = current->value;
    current = current->next;
  }
  return unpacked;
}

// ----------------------------------------------------
// Mathematical Operations
// ----------------------------------------------------

CobMatrix *CobMatrix_compose(CobMatrix *this_m, CobMatrix *that_m) {
  // Assert target of that_m == source of this_m
  // that_m applies first: that_m->source -> that_m->target
  // this_m applies second: this_m->source -> this_m->target
  // In Java "this.compose(matrix)" implies "this * matrix".

  CobMatrix *result = CobMatrix_create(that_m->source, this_m->target, false);

  for (int i = 0; i < this_m->target->n; i++) {
    MatrixEntry *rowI = this_m->entries[i].head;
    while (rowI != NULL) {
      int j = rowI->column_index;
      if (j < 0 || j >= that_m->target->n) {
        kh_fatal("invalid row index in matrix composition");
      }
      // Iterate over that_m's row j
      MatrixEntry *that_rowJ = that_m->entries[j].head;
      while (that_rowJ != NULL) {
        int k = that_rowJ->column_index;
        if (k < 0 || k >= result->source->n) {
          kh_fatal("invalid column index in matrix composition");
        }
        LCCC *composed = LCCC_compose(rowI->value, that_rowJ->value);
        if (composed != NULL) {
          // addEntry consumes composed, including the zero case.
          CobMatrix_addEntry(result, i, k, composed);
        }
        that_rowJ = that_rowJ->next;
      }
      rowI = rowI->next;
    }
  }
  return result;
}

void CobMatrix_multiply(CobMatrix *m, RingElement *n) {
  // Multiply all entries by n
  for (int i = 0; i < m->target->n; i++) {
    MatrixEntry *current = m->entries[i].head;
    while (current != NULL) {
      LCCC *mult = LCCC_multiply(current->value, n);
      // Replace value
      LCCC_free(current->value);
      current->value = mult;
      current = current->next;
    }
  }
}

void CobMatrix_add(CobMatrix *dest, CobMatrix *src) {
  // Both must have same source and target n.
  assert(dest->target->n == src->target->n);
  for (int i = 0; i < src->target->n; i++) {
    MatrixEntry *src_entry = src->entries[i].head;
    while (src_entry != NULL) {
      // addEntry takes ownership, so clone borrowed values from src.
      CobMatrix_addEntry(dest, i, src_entry->column_index,
                         LCCC_clone(src_entry->value));
      src_entry = src_entry->next;
    }
  }
}

void CobMatrix_reduce(CobMatrix *m) {
  for (int i = 0; i < m->target->n; i++) {
    MatrixRow *row = &m->entries[i];
    MatrixEntry *prev = NULL;
    MatrixEntry *current = row->head;

    while (current != NULL) {
      LCCC *old_value = current->value;
      LCCC *reduced = LCCC_reduce(old_value);
      LCCC_free(old_value);

      if (reduced == NULL || LCCC_isZero(reduced)) {
        MatrixEntry *to_delete = current;
        if (prev) {
          prev->next = current->next;
        } else {
          row->head = current->next;
        }
        current = current->next;
        LCCC_free(reduced);
        free(to_delete);
      } else {
        current->value = reduced;
        prev = current;
        current = current->next;
      }
    }
  }
}

bool CobMatrix_isZero(CobMatrix *m) {
  for (int i = 0; i < m->target->n; i++) {
    MatrixEntry *current = m->entries[i].head;
    while (current != NULL) {
      if (current->value != NULL && !LCCC_isZero(current->value)) {
        return false;
      }
      current = current->next;
    }
  }
  return true;
}

// ----------------------------------------------------
// Row / Column Extraction
// ----------------------------------------------------
/*
 * Extract one source column from m.
 *
 * This mutates m in place by removing the selected column and reindexing
 * remaining entries. The returned matrix is a one-column slice using shallow
 * SmoothingColumn data and moved LCCC pointers.
 */
CobMatrix *CobMatrix_extractColumn(CobMatrix *m, int column_idx) {
  SmoothingColumn *newTarget = m->target;
  SmoothingColumn *newSource =
      SmoothingColumn_create();
  newSource->n = 1;
  newSource->numbers = (int *)kh_malloc(sizeof(int));
  newSource->smoothings = (Cap **)kh_malloc(sizeof(Cap *));

  // Copy the removed column's smoothing column info
  newSource->numbers[0] = m->source->numbers[column_idx];
  newSource->smoothings[0] = Cap_retain(m->source->smoothings[column_idx]);

  CobMatrix *res = CobMatrix_create(newSource, newTarget, true);
  SmoothingColumn_free(newSource);

  for (int i = 0; i < m->target->n; i++) {
    MatrixEntry *prev = NULL;
    MatrixEntry *current = m->entries[i].head;

    while (current != NULL) {
      if (current->column_index == column_idx) {
        CobMatrix_putEntry(res, i, 0, current->value);

        // Remove from m
        if (prev)
          prev->next = current->next;
        else
          m->entries[i].head = current->next;

        MatrixEntry *next = current->next;
        free(current);
        current = next;
        continue;
      } else if (current->column_index > column_idx) {
        // decrement indices to shift left
        current->column_index--;
      }
      prev = current;
      current = current->next;
    }
  }

  // Create a new source column to avoid mutating shared arrays
  SmoothingColumn *new_m_source =
      SmoothingColumn_create();
  new_m_source->n = m->source->n - 1;
  if (new_m_source->n > 0) {
    new_m_source->numbers = (int *)kh_malloc(sizeof(int) * new_m_source->n);
    new_m_source->smoothings = (Cap **)kh_malloc(sizeof(Cap *) * new_m_source->n);
    for (int i = 0, j = 0; i < m->source->n; i++) {
      if (i == column_idx)
        continue;
      new_m_source->numbers[j] = m->source->numbers[i];
      new_m_source->smoothings[j] = Cap_retain(m->source->smoothings[i]);
      j++;
    }
  } else {
    new_m_source->numbers = NULL;
    new_m_source->smoothings = NULL;
  }

  SmoothingColumn_free(m->source);
  m->source = new_m_source;

  return res;
}
/*
 * Extract one target row from m.
 *
 * This mutates m in place by removing the selected row and shifting remaining
 * rows upward. The returned matrix is a one-row slice using shallow
 * SmoothingColumn data and moved LCCC pointers.
 */
CobMatrix *CobMatrix_extractRow(CobMatrix *m, int row_idx) {
  SmoothingColumn *newSource = m->source;
  SmoothingColumn *newTarget =
      SmoothingColumn_create();
  newTarget->n = 1;
  newTarget->numbers = (int *)kh_malloc(sizeof(int));
  newTarget->smoothings = (Cap **)kh_malloc(sizeof(Cap *));

  // Copy the removed row's smoothing target info
  newTarget->numbers[0] = m->target->numbers[row_idx];
  newTarget->smoothings[0] = Cap_retain(m->target->smoothings[row_idx]);

  CobMatrix *res = CobMatrix_create(newSource, newTarget, true);
  SmoothingColumn_free(newTarget);

  // Move list out
  res->entries[0].head = m->entries[row_idx].head;
  m->entries[row_idx].head = NULL;

  // Shift remaining rows up
  for (int i = row_idx + 1; i < m->target->n; i++) {
    m->entries[i - 1] = m->entries[i];
  }
  m->entries[m->target->n - 1].head = NULL;

  // Create a new target column to avoid mutating shared arrays
  SmoothingColumn *new_m_target =
      SmoothingColumn_create();
  new_m_target->n = m->target->n - 1;
  if (new_m_target->n > 0) {
    new_m_target->numbers = (int *)kh_malloc(sizeof(int) * new_m_target->n);
    new_m_target->smoothings = (Cap **)kh_malloc(sizeof(Cap *) * new_m_target->n);
    for (int i = 0, j = 0; i < m->target->n; i++) {
      if (i == row_idx)
        continue;
      new_m_target->numbers[j] = m->target->numbers[i];
      new_m_target->smoothings[j] = Cap_retain(m->target->smoothings[i]);
      j++;
    }
  } else {
    new_m_target->numbers = NULL;
    new_m_target->smoothings = NULL;
  }

  SmoothingColumn_free(m->target);
  m->target = new_m_target;

  return res;
}

SmoothingColumn *SmoothingColumn_clone(SmoothingColumn *col) {
  if (col == NULL) return NULL;
  SmoothingColumn *clone = SmoothingColumn_create();

  clone->n = col->n;
  if (col->n <= 0) {
    clone->numbers = NULL;
    clone->smoothings = NULL;
    return clone;
  }
  clone->numbers = (int *)kh_malloc((size_t)col->n * sizeof(int));
  clone->smoothings = (Cap **)kh_malloc((size_t)col->n * sizeof(Cap *));

  memcpy(clone->numbers, col->numbers, (size_t)col->n * sizeof(int));
  for (int i = 0; i < col->n; ++i)
    clone->smoothings[i] = Cap_retain(col->smoothings[i]);
  return clone;
}
bool SmoothingColumn_equals(SmoothingColumn *a, SmoothingColumn *b) {
  if (a == b) return true;
  if (a == NULL || b == NULL || a->n != b->n) return false;
  for (int i = 0; i < a->n; i++) {
    if (a->numbers[i] != b->numbers[i]) return false;
    if (!Cap_equals(a->smoothings[i], b->smoothings[i])) return false;
  }
  return true;
}
