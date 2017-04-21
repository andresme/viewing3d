#ifndef POLYGONS_VERTEX_H
#define POLYGONS_VERTEX_H

#include "../util/math/mathUtils.h"

typedef struct vertex {
    long double x;
    long double y;
    long double z;
} vertex;

typedef struct polygon {
    int vertices[3];
} polygon;

matrix applyTransformation(vertex v, matrix transformation);
vertex applyAll(vertex v, matrix mper, matrix mvv3dv);

#endif //POLYGONS_VERTEX_H
