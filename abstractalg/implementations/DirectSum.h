#ifndef DIRECTSUM_H
#define DIRECTSUM_H

#include <stdbool.h>

//Purpose: A container representing a mathematical direct sum
//Method: Wraps a dynamic array of pointers to algebraic objects
typedef struct {
    void** elements;
    int size;
    int capacity;
} DirectSum;

//Purpose: Initialize a DirectSum container
//Returns: true if memory allocation succeeded, false otherwise
bool DirectSum_init(DirectSum* ds);

//Purpose: Add an object to the direct sum
//Returns: true if successful, false if memory reallocation failed
bool DirectSum_add(DirectSum* ds, void* element);

//Purpose: Frees the dynamically allocated array to prevent memory leaks
void DirectSum_free(DirectSum* ds);

#endif