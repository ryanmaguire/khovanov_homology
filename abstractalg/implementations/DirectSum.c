#include "DirectSum.h"
#include <stdlib.h>

//Purpose: Initialize a DirectSum container
//Arguments: ds
//Argument descriptions: ds is a pointer to the DirectSum to initialize
//Output: void
//Method: Allocates initial memory for the elements array and sets size to 0
void DirectSum_init(DirectSum* ds) {
    ds->size = 0;
    ds->capacity = 10;
    ds->elements = malloc(sizeof(void*) * ds->capacity);
}

//Purpose: Add an object to the direct sum
//Arguments: ds, element
//Argument descriptions: ds is the container, element is the pointer to the object
//Output: void
//Method: Checks capacity, resizes if necessary, and appends the element pointer
void DirectSum_add(DirectSum* ds, void* element) {
    if (ds->size >= ds->capacity) {
        ds->capacity *= 2;
        ds->elements = realloc(ds->elements, sizeof(void*) * ds->capacity);
    }
    ds->elements[ds->size++] = element;
}