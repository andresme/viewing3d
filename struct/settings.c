#include <stdlib.h>
#include "settings.h"



long double* createVector(long double vector[], int size) {
    long double *assignedVector = malloc(size * sizeof(long double));
    for(int i = 0; i < size; i++) {
        assignedVector[i] = vector[i];
    }
    return assignedVector;
}