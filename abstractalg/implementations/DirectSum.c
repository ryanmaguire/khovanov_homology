#include "DirectSum.h"
#include <stdlib.h>

//Purpose: Initialize a DirectSum container
//Arguments: ds - pointer to the DirectSum to initialize
//Returns: true if memory allocation succeeded, false otherwise
bool DirectSum_init(DirectSum* ds) {
    //safety check
    if (!ds) return false; 
    
    ds->size = 0;
    ds->capacity = 10;
    ds->elements = malloc(sizeof(void*) * ds->capacity);
    
    //return true if malloc succeeded, false if we ran out of RAM
    return ds->elements != NULL; 
}

//Purpose: Add an object to the direct sum
//Arguments: ds - the container, element - the pointer to the object
//Returns: true if successful, false if memory reallocation failed
bool DirectSum_add(DirectSum* ds, void* element) {
    if (!ds) return false;
    
    if (ds->size >= ds->capacity) {
        int new_capacity = ds->capacity > 0 ? ds->capacity * 2 : 10;
        void** new_elements = realloc(ds->elements, sizeof(void*) * new_capacity);
        
        //safety
        if (!new_elements) {
            return false; 
        }
        
        ds->elements = new_elements;
        ds->capacity = new_capacity;
    }
    
    ds->elements[ds->size++] = element;
    return true;
}

//Purpose: Frees the dynamically allocated array to prevent memory leaks
//Arguments: ds - pointer to the DirectSum to clean up
//Returns: void
void DirectSum_free(DirectSum* ds) {
    if (!ds) return;
    
    if (ds->elements) {
        free(ds->elements);
        ds->elements = NULL; 
    }
    
    ds->size = 0;
    ds->capacity = 0;
}