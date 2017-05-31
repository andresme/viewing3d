#include "mathUtils.h"
#include <math.h>
#include <stdlib.h>

matrix initMatrix(int height, int width) {
    matrix m;
    m.height = height;
    m.width = width;
    m.values = malloc(m.height * sizeof(long double*));
    for(int i = 0; i < m.height; i++) {
        m.values[i] = malloc(m.width * sizeof (long double));
    }
    for(int i = 0; i < m.height; i++) {
        for(int j = 0; j < m.width; j++) {
            m.values[i][j] = 0;
        }
    }
    return m;
}

matrix initVector(int size) {
    return initMatrix(1, size);
}

matrix identityMatrix(int size) {
    matrix result = initMatrix(size, size);

    for(int i = 0; i < result.width; i++) {
            result.values[i][i] = 1;
    }
    return result;
}

matrix multiplyMatrixByMatrix(matrix m1, matrix m2) {
    matrix result = initMatrix(m1.height, m2.width);

    for(int i = 0; i < m1.height; i++) {
        for(int j = 0; j < m2.width; j++) {
            for(int k = 0; k < m1.width; k++) {
                result.values[i][j] += m1.values[i][k] * m2.values[k][j];
            }
        }
    }
    return result;
}

long double getVectorLength(matrix v) {
    long double sum = 0;
    for(int i = 0; i < v.width; i++) {
        sum += v.values[0][i] * v.values[0][i];
    }
    return sqrt(sum);
}

matrix normalizeVector(matrix v) {
    matrix result = initMatrix(v.height, v.width);
    long double norm = getVectorLength(v);
    for(int i = 0; i < v.width; i++) {
        result.values[0][i] = v.values[0][i] / norm;
    }
    return result;
}

matrix crossProduct(matrix u, matrix v) {
    matrix result = initMatrix(1, 3);
    result.values[0][0] = u.values[0][1] * v.values[0][2] - u.values[0][2] * v.values[0][1];
    result.values[0][1] = u.values[0][2] * v.values[0][0] - u.values[0][0] * v.values[0][2];
    result.values[0][2] = u.values[0][0] * v.values[0][1] - u.values[0][1] * v.values[0][0];

    return result;
}

matrix getHomogeneousVector(matrix v1) {
    matrix result = initMatrix(1, 4);

    for(int i = 0; i < 3; i++){
        result.values[0][i] = v1.values[0][i];
    }
    result.values[0][3] = 1;
    return result;

}

matrix translateOriginMatrix(matrix v) {
    matrix result = identityMatrix(4);

    for(int i = 0; i < 3; i++) {
        result.values[i][3] = - v.values[0][i];
    }

    return result;
}

matrix translate(matrix v) {
    matrix result = identityMatrix(4);

    for(int i = 0; i < 3; i++) {
        result.values[i][3] = v.values[0][i];
    }

    return result;
}

matrix substract(matrix v1, matrix v2) {
    matrix result = initMatrix(1, v1.width);
    for(int i = 0; i < v1.width; i++) {
        result.values[0][i] = v1.values[0][i] - v2.values[0][i];
    }
    return result;
}

matrix dopShearMatrix(matrix v) {
    matrix result = identityMatrix(4);
    long double shx = - v.values[0][0] / v.values[0][2];
    long double shy = - v.values[0][1] / v.values[0][2];

    result.values[0][2] = shx;
    result.values[1][2] = shy;

    return result;

}

matrix scaleMatrix(matrix v) {
    matrix result = identityMatrix(4);
    for(int i  = 0; i < 3; i++) {
        result.values[i][i] = v.values[0][i];
    }
    return result;
}

matrix transpose(matrix m) {
    matrix result = initMatrix(m.width, m.height);

    for(int i = 0; i < m.width; i++) {
        for(int j = 0; j < m.height; j++) {
            result.values[i][j] = m.values[j][i];
        }
    }
    return result;
}

void freeMatrix(matrix m) {
    for(int i = 0; i < m.height; i++){
        free(m.values[i]);
    }
    free(m.values);
}