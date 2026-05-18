#ifndef MORPHISM_H
#define MORPHISM_H

#include <stdbool.h>

typedef struct Morphism Morphism;

//compares two morphisms for equality
bool Morphism_equals(Morphism* m1, Morphism* m2);

#endif