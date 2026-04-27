#ifndef MORPHISM_H
#define MORPHISM_H

typedef struct Morphism Morphism;

struct Morphism {
    //Purpose: Get a string representation of the morphism
    //Arguments: self
    //Argument descriptions: self is a pointer to the Morphism instance
    //Output: char*
    //Output description:an allocated string
    char* (*toString)(Morphism* self);

    //Purpose: Compare two morphisms for equality
    //Arguments: self, other
    //Argument descriptions: pointers to the two morphisms to compare
    //Output: int
    //Output description: 1 if equal, 0 otherwise
    int (*equals)(Morphism* self, Morphism* other);

    //Purpose: Compose this morphism with another
    //Arguments: self, other
    //Argument descriptions: self and other are the morphisms being composed
    //Output: Morphism*
    //Output description: A new morphism representing the composition
    Morphism* (*compose)(Morphism* self, Morphism* other);

    int (*compareTo)(struct Morphism* self, struct Morphism* other);
};

#endif
