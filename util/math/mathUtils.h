#ifndef POLYGONS_MATHUTILS_H
#define POLYGONS_MATHUTILS_H

typedef struct matrix {
    int width;
    int height;
    long double **values;
} matrix;

matrix initMatrix(int height, int width);
matrix initVector(int size);
matrix identityMatrix(int size);
matrix multiplyMatrixByMatrix(matrix m1, matrix m2);
long double getVectorLength(matrix v);
matrix normalizeVector(matrix v);
matrix crossProduct(matrix u, matrix v);
matrix getHomogeneousVector(matrix v1);
matrix translateOriginMatrix(matrix v);
matrix translate(matrix v);
matrix substract(matrix v1, matrix v2);
matrix dopShearMatrix(matrix v);
matrix scaleMatrix(matrix v);
void freeMatrix(matrix m);
matrix transpose(matrix m);

#endif //POLYGONS_MATHUTILS_H
