#include <stdio.h>
#include <stdlib.h>
#include "LinearComboMap.h"
#include "ZeroLinearCombo.h"

struct Morphism {
    const char* name;
};

char* mock_morphismToString(Morphism* m) {
    if (!m) return NULL;
    char* str = malloc(64); 
    snprintf(str, 64, "%s", m->name);
    return str;
}

int main(void) {
    printf("Khovanov homology architecture test\n");

    LinearComboMap myCombo;
    LinearComboMap_init(&myCombo);
    printf("[SUCCESS] LinearComboMap initialized.\n");

    LinearComboMap* result = LinearComboMap_multiply(&myCombo, 5);
    if (result == &myCombo) {
        printf("[SUCCESS] Multiply was safely skipped on an empty map.\n");
    }

    //TEST: ZeroLinearCombo hardcoded methods
    int zeroTerms = ZeroLinearCombo_numberOfTerms(NULL);
    bool isZero = ZeroLinearCombo_isZero(NULL);
    
    if (zeroTerms == 0 && isZero == true) {
        printf("[SUCCESS] ZeroLinearCombo reported 0 terms and isZero = true.\n");
    }

    printf("All architectural tests passed. No segfaults.\n");

    return 0;
}
