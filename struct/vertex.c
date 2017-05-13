#include <malloc.h>
#include "vertex.h"
#include "../util/math/mathUtils.h"

matrix vertexToHomogeneousVector(vertex v) {
    matrix result = initVector(4);
    result.values[0][0] = v.x;
    result.values[0][1] = v.y;
    result.values[0][2] = v.z;
    result.values[0][3] = 1;
    return result;
}

matrix applyTransformation(vertex v, matrix transformation) {
    matrix homogeneous = vertexToHomogeneousVector(v);
    matrix result = multiplyMatrixByMatrix(transformation, transpose(homogeneous));
    return transpose(result);
}

vertex applyAll(vertex v, matrix mper, matrix mvv3dv) {
    matrix transformed = transpose(applyTransformation(v, identityMatrix(4)));
    transformed = multiplyMatrixByMatrix(mper, transformed);

    transformed = multiplyMatrixByMatrix(mvv3dv, transformed);
    v.x = transformed.values[0][0] / transformed.values[3][0];
    v.y = transformed.values[1][0] / transformed.values[3][0];
    v.z = -transformed.values[3][0] * mvv3dv.values[2][3];
    return v;
}