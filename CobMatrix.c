#include "CobMatrix.h"
#include <assert.h>
#include <stdlib.h>

// Helper to create a new linked list node
static MatrixEntry *create_entry(int col_idx, LCCC *value) {
  MatrixEntry *entry = (MatrixEntry *)malloc(sizeof(MatrixEntry));
  entry->column_index = col_idx;
  entry->value = value;
  entry->next = NULL;
  return entry;
}

// Helper to free a matrix row list (the row itself is a value struct usually)
static void free_matrix_row(MatrixRow *row) {
  MatrixEntry *current = row->head;
  while (current != NULL) {
    MatrixEntry *next = current->next;
    // Optionally free the LCCC here if ownership dictates:
    // LCCC_free(current->value);
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
  CobMatrix *m = (CobMatrix *)malloc(sizeof(CobMatrix));
  if (shared) {
    m->source = source;
    m->target = target;
  } else {
    m->source = SmoothingColumn_clone(source);
    m->target = SmoothingColumn_clone(target);
  }

  m->entries = (MatrixRow *)calloc(m->target->n, sizeof(MatrixRow));
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
  // Note: ownership of source/target is murky here depending on 'shared',
  // assume handled elsewhere or add logic.
  free(m);
}

// ----------------------------------------------------
// Core Matrix Operations
// ----------------------------------------------------

void CobMatrix_putEntry(CobMatrix *m, int row_idx, int col_idx, LCCC *lc) {
  if (lc == NULL || LCCC_isZero(lc)) {
    return;
  }

  MatrixRow *row = &m->entries[row_idx];
  MatrixEntry *new_entry = create_entry(col_idx, lc);

  // Insert into linked list (keep sorted by col_idx ideally, but appending is
  // fine for now) Here we append at head for simplicity, but for matrix
  // multiplications sorted order is better.
  new_entry->next = row->head;
  row->head = new_entry;
}

void CobMatrix_addEntry(CobMatrix *m, int row_idx, int col_idx, LCCC *t) {
  if (t == NULL || LCCC_isZero(t)) {
    return;
  }

  MatrixRow *row = &m->entries[row_idx];
  MatrixEntry *prev = NULL;
  MatrixEntry *current = row->head;

  while (current != NULL) {
    if (current->column_index == col_idx) {
      LCCC *sum = LCCC_add(current->value, t);
      if (LCCC_isZero(sum)) {
        // remove entry
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

  // Not found, append
  MatrixEntry *new_entry = create_entry(col_idx, t);
  new_entry->next = row->head;
  row->head = new_entry;
}

LCCC **CobMatrix_unpackRow(CobMatrix *m, int row_idx) {
  LCCC **unpacked = (LCCC **)calloc(m->source->n, sizeof(LCCC *));
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
      // Iterate over that_m's row j
      MatrixEntry *that_rowJ = that_m->entries[j].head;
      while (that_rowJ != NULL) {
        int k = that_rowJ->column_index;
        LCCC *composed = LCCC_compose(rowI->value, that_rowJ->value);

        if (composed != NULL && !LCCC_isZero(composed)) {
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
      CobMatrix_addEntry(dest, i, src_entry->column_index, src_entry->value);
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
      LCCC *reduced = LCCC_reduce(current->value);
      if (reduced == NULL || LCCC_isZero(reduced)) {
        // remove entry
        MatrixEntry *to_delete = current;
        if (prev) {
          prev->next = current->next;
        } else {
          row->head = current->next;
        }
        current = current->next;
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

CobMatrix *CobMatrix_extractColumn(CobMatrix *m, int column_idx) {
  SmoothingColumn *newTarget = m->target;
  SmoothingColumn *newSource =
      (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  newSource->n = 1;
  newSource->numbers = (int *)malloc(sizeof(int));
  newSource->smoothings = (Cap **)malloc(sizeof(Cap *));

  // Copy the removed column's smoothing column info
  newSource->numbers[0] = m->source->numbers[column_idx];
  newSource->smoothings[0] = m->source->smoothings[column_idx];

  CobMatrix *res = CobMatrix_create(newSource, newTarget, true);

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
      (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  new_m_source->n = m->source->n - 1;
  if (new_m_source->n > 0) {
    new_m_source->numbers = (int *)malloc(sizeof(int) * new_m_source->n);
    new_m_source->smoothings = (Cap **)malloc(sizeof(Cap *) * new_m_source->n);
    for (int i = 0, j = 0; i < m->source->n; i++) {
      if (i == column_idx)
        continue;
      new_m_source->numbers[j] = m->source->numbers[i];
      new_m_source->smoothings[j] = m->source->smoothings[i];
      j++;
    }
  } else {
    new_m_source->numbers = NULL;
    new_m_source->smoothings = NULL;
  }

  // Note: we intentionally do not free the old m->source array
  // because it may be shared if shared=true.
  m->source = new_m_source;

  return res;
}

CobMatrix *CobMatrix_extractRow(CobMatrix *m, int row_idx) {
  SmoothingColumn *newSource = m->source;
  SmoothingColumn *newTarget =
      (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  newTarget->n = 1;
  newTarget->numbers = (int *)malloc(sizeof(int));
  newTarget->smoothings = (Cap **)malloc(sizeof(Cap *));

  // Copy the removed row's smoothing target info
  newTarget->numbers[0] = m->target->numbers[row_idx];
  newTarget->smoothings[0] = m->target->smoothings[row_idx];

  CobMatrix *res = CobMatrix_create(newSource, newTarget, true);

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
      (SmoothingColumn *)malloc(sizeof(SmoothingColumn));
  new_m_target->n = m->target->n - 1;
  if (new_m_target->n > 0) {
    new_m_target->numbers = (int *)malloc(sizeof(int) * new_m_target->n);
    new_m_target->smoothings = (Cap **)malloc(sizeof(Cap *) * new_m_target->n);
    for (int i = 0, j = 0; i < m->target->n; i++) {
      if (i == row_idx)
        continue;
      new_m_target->numbers[j] = m->target->numbers[i];
      new_m_target->smoothings[j] = m->target->smoothings[i];
      j++;
    }
  } else {
    new_m_target->numbers = NULL;
    new_m_target->smoothings = NULL;
  }

  // Note: we intentionally do not free the old m->target array
  // because it may be shared if shared=true.
  m->target = new_m_target;

  return res;
}
