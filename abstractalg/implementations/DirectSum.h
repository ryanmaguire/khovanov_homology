#ifndef DIRECTSUM_H
#define DIRECTSUM_H

//Purpose: A container representing a mathematical direct sum
//Method: Wraps a dynamic array of pointers to algebraic objects
typedef struct {
    void** elements;
    int size;
    int capacity;
} DirectSum;

//Purpose: Initialize a DirectSum container
//Arguments: ds
//Argument descriptions: ds is a pointer to the DirectSum to initialize
//Output type: void
void DirectSum_init(DirectSum* ds);

//Purpose: Add an object to the direct sum
//Arguments: ds, element
//Argument descriptions: ds is the container, element is the pointer to the object
//Output type: void
void DirectSum_add(DirectSum* ds, void* element);

#endif