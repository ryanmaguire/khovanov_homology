#ifndef RING_H
#define RING_H

typedef struct Ring Ring;
struct Ring {
    //Purpose: Add two elements of the ring
    //Arguments: a, b
    //Argument descriptions: pointers to the ring elements to add
    //Output: void*
    //Output description: A pointer to the resulting ring element
    void* (*add)(void* a, void* b);

    //Purpose: Multiply two elements of the ring
    //Arguments: a, b
    //Argument descriptions: pointers to the ring elements to multiply
    //Output: void*
    //Output description: A pointer to the resulting ring element
    void* (*multiply)(void* a, void* b);

    //Purpose: check if the element is the zero of the ring
    //Arguments: element
    //Argument descriptions: the element to check
    //Output: int
    //Output description: 1 if zero, 0 otherwise
    int (*isZero)(void* element);
};

#endif