#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <string.h>
#include "util/math/mathUtils.h"
#include "struct/vertex.h"
#include "struct/settings.h"
#include <math.h>
#include <values.h>


long double WHITE[] = {1, 1, 1};
long double BLACK[] = {0, 0, 0};

void printMatrix(matrix m, const char *name) {
    printf("======%s======\n", name);
    for(int i = 0; i < m.height; i++) {
        for(int j = 0; j < m.width; j++) {
            printf("%Lf\t",m.values[i][j]);
        }
        printf("\n");
    }
}

int verticesCount;
int facesCount;
int type = 0;

long double frameBuffer[500][500][3];
long double zbuffer[500][500];

vertex *vertices;
vertex *transformedVertices;
polygon *faces;

scene_settings settings;

void bres(int x0, int y0, long double z0, int x1, int y1, long double z1, long double normal1[], long double normal2[]);
long double interpolate(long double x, long double x0, long double x1, long double y0, long double y1);

void clearBuffer() {
    for(int i = 0; i < 500; i++) {
        for(int j = 0; j < 500; j++) {
            frameBuffer[i][j][0] = WHITE[0];
            frameBuffer[i][j][1] = WHITE[1];
            frameBuffer[i][j][2] = WHITE[2];
            zbuffer[i][j] = -INFINITY;
        }
    }
}

void draw_pixel(int x, int y, long double z, long double color[]) {
    if(x >= 0 && x <= 500 && y >= 0 && y <= 500) {
        if(z > zbuffer[x][y]){
            frameBuffer[x][y][0] = color[0];
            frameBuffer[x][y][1] = color[1];
            frameBuffer[x][y][2] = color[2];
            zbuffer[x][y] = z;
        }

    }
}

vertex* getOrderedVertices(vertex vertex1, vertex vertex2, vertex vertex3) {
    vertex* result = malloc(3 * sizeof(vertex));

    if(vertex1.y < vertex2.y) {
        if(vertex1.y < vertex3.y) {
            result[0] = vertex1;
            if(vertex2.y < vertex3.y) {
                result[1] = vertex2;
                result[2] = vertex3;
            } else {
                result[1] = vertex3;
                result[2] = vertex2;
            }
        } else {
            result[0] = vertex3;
            result[1] = vertex1;
            result[2] = vertex2;
        }
    } else if(vertex2.y < vertex3.y) {
        result[0] = vertex2;
        if(vertex1.y < vertex3.y) {
            result[1] = vertex1;
            result[2] = vertex3;
        } else {
            result[1] = vertex3;
            result[2] = vertex1;
        }
    } else {
        result[0] = vertex3;
        result[1] = vertex2;
        result[2] = vertex1;
    }

    return result;
}

void drawBottom(vertex v1, vertex v2, vertex v3) {

    long double invslope1 = (v2.x - v1.x) / (v2.y - v1.y);
    long double invslope2 = (v3.x - v1.x) / (v3.y - v1.y);

    long double curx1 = v1.x;
    long double curx2 = v1.x;

    long double curz1 = v1.z;
    long double curz2 = v1.z;

    long double z_inc_left = (v2.z - v1.z)/(v2.y - v1.y);
    long double z_inc_right = (v3.z - v1.z)/(v3.y - v1.y);

    for (long double scanlineY = v1.y; scanlineY <= v2.y; scanlineY++)
    {
        long double normal1[3] = {
            interpolate(scanlineY, v1.y, v2.y, v1.avg_normal[0], v2.avg_normal[0]),
            interpolate(scanlineY, v1.y, v2.y, v1.avg_normal[1], v2.avg_normal[1]),
            interpolate(scanlineY, v1.y, v2.y, v1.avg_normal[2], v2.avg_normal[2]),
        };

        long double normal2[3] = {
                interpolate(scanlineY, v1.y, v3.y, v1.avg_normal[0], v3.avg_normal[0]),
                interpolate(scanlineY, v1.y, v3.y, v1.avg_normal[1], v3.avg_normal[1]),
                interpolate(scanlineY, v1.y, v3.y, v1.avg_normal[2], v3.avg_normal[2]),
        };

        bres(curx1, scanlineY, curz1, curx2, scanlineY, curz2, normal1, normal2);
        curx1 += invslope1;
        curx2 += invslope2;

        curz1 += z_inc_left;
        curz2 += z_inc_right;
    }

    bres(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v1.avg_normal, v2.avg_normal);
    bres(v2.x, v2.y, v2.z, v3.x, v3.y, v3.z, v2.avg_normal, v3.avg_normal);
    bres(v3.x, v3.y, v3.z, v1.x, v1.y, v1.z, v3.avg_normal, v1.avg_normal);
}

void drawTop(vertex v1, vertex v2, vertex v3) {
    long double invslope1 = (v3.x - v1.x) / (v3.y - v1.y);
    long double invslope2 = (v3.x - v2.x) / (v3.y - v2.y);

    long double curx1 = v3.x;
    long double curx2 = v3.x;

    long double curz1 = v3.z;
    long double curz2 = v3.z;

    long double z_inc_left = (v1.z - v3.z)/(v1.y - v3.y);
    long double z_inc_right = (v2.z - v3.z)/(v2.y - v3.y);

    for (long double scanlineY = v3.y; scanlineY >= v1.y; scanlineY--) {
        long double normal1[3] = {
                interpolate(scanlineY, v3.y, v1.y, v3.avg_normal[0], v1.avg_normal[0]),
                interpolate(scanlineY, v3.y, v1.y, v3.avg_normal[1], v1.avg_normal[1]),
                interpolate(scanlineY, v3.y, v1.y, v3.avg_normal[2], v1.avg_normal[2]),
        };

        long double normal2[3] = {
                interpolate(scanlineY, v3.y, v2.y, v3.avg_normal[0], v2.avg_normal[0]),
                interpolate(scanlineY, v3.y, v2.y, v3.avg_normal[1], v2.avg_normal[1]),
                interpolate(scanlineY, v3.y, v2.y, v3.avg_normal[2], v2.avg_normal[2]),
        };
        bres(curx1, scanlineY, curz1, curx2, scanlineY, curz2, normal1, normal2);
        curx1 -= invslope1;
        curx2 -= invslope2;

        curz1 += z_inc_left;
        curz2 += z_inc_right;
    }
    bres(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v1.avg_normal, v2.avg_normal);
    bres(v2.x, v2.y, v2.z, v3.x, v3.y, v3.z, v2.avg_normal, v3.avg_normal);
    bres(v3.x, v3.y, v3.z, v1.x, v1.y, v1.z, v3.avg_normal, v1.avg_normal);
}


long double* getNormal(vertex v1, vertex v2, vertex v3) {
    long double *result = malloc(3 * sizeof(long double));
    /**
    * V = o[1] - o[0]
    * W = o[2] - o[0]
    * Nx=(Vy∗Wz)−(Vz∗Wy)
      Ny=(Vz∗Wx)−(Vx∗Wz)
      Nz=(Vx∗Wy)−(Vy∗Wx)
    */
    long double V[3] = {
            v2.x - v1.x,
            v2.y - v1.y,
            v2.z - v1.z
    };

    long double W[3] = {
            v3.x - v1.x,
            v3.y - v1.y,
            v3.z - v1.z
    };
    result[0] = V[1]*W[2] - V[2]*W[1];
    result[1] = V[2]*W[0] - V[0]*W[2];
    result[2] = V[0]*W[1] - V[1]*W[0];

    long double size = sqrt(pow((double) result[0], 2) +
                            pow((double) result[1], 2) +
                            pow((double) result[2], 2));

    result[0] = result[0] / size;
    result[1] = result[1] / size;
    result[2] = result[2] / size;

    return result;
}

void drawTriangle(vertex vertex1, vertex vertex2, vertex vertex3) {
    vertex* ordered = getOrderedVertices(vertex1, vertex2, vertex3);

    if (ordered[1].y == ordered[2].y) {
        drawBottom(ordered[0], ordered[1], ordered[2]);
    } else if (ordered[0].y == ordered[1].y) {
        drawTop(ordered[0], ordered[1], ordered[2]);
    } else {
        vertex vertex4 = {(ordered[0].x + ((ordered[1].y - ordered[0].y) / (ordered[2].y - ordered[0].y)) * (ordered[2].x - ordered[0].x)),
                          ordered[1].y,
                          (ordered[0].z * (ordered[2].y - ordered[1].y) + ordered[2].z * (ordered[1].y - ordered[0].y)) / (ordered[2].y - ordered[0].y)};

        vertex4.avg_normal[0] = interpolate(vertex4.y, ordered[0].y, ordered[2].y, ordered[0].avg_normal[0], ordered[2].avg_normal[0]);
        vertex4.avg_normal[1] = interpolate(vertex4.y, ordered[0].y, ordered[2].y, ordered[0].avg_normal[1], ordered[2].avg_normal[1]);
        vertex4.avg_normal[2] = interpolate(vertex4.y, ordered[0].y, ordered[2].y, ordered[0].avg_normal[2], ordered[2].avg_normal[2]);

        drawBottom(ordered[0], ordered[1], vertex4);
        drawTop(ordered[1], vertex4, ordered[2]);
    }

}

long double interpolate(long double x, long double x0, long double x1, long double y0, long double y1) {
    return (y0*(x1-x)+y1*(x-x0))/(x1-x0);
}

void bres(int x0, int y0, long double z0, int x1, int y1, long double z1, long double normal1[], long double normal2[]) {
    int original = x0;
    int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
    int dy = abs(y1-y0), sy = y0<y1 ? 1 : -1;
    int err = (dx>dy ? dx : -dy)/2, e2;

    for(;;){
        long double current_z = interpolate(x0, original, x1, z0, z1);
        long double light_vector[3] = {250 - x0, 250 - y0, 100 - current_z};

        long double size = sqrt(pow(light_vector[0],2) + pow(light_vector[1], 2) + pow(light_vector[2],2));

        light_vector[0] = light_vector[0] / size;
        light_vector[1] = light_vector[1] / size;
        light_vector[2] = light_vector[2] / size;

        long double pixel_normal[3] = {
                interpolate(x0, original, x1, normal1[0], normal2[0]),
                interpolate(x0, original, x1, normal1[1], normal2[1]),
                interpolate(x0, original, x1, normal1[2], normal2[2]),
        };

        long double dot_product =
                pixel_normal[0] * light_vector[0] +
                pixel_normal[1] * light_vector[1] +
                pixel_normal[2] * light_vector[2];

        long double intensity = sqrt(pow(dot_product, 2));

        long double color[3] = {
                intensity*114.0/255.0,
                intensity*64.0/255.0,
                intensity*11.0/255.0,
        };
        draw_pixel(x0,y0, current_z, color);
        if (x0==x1 && y0==y1) break;
        e2 = err;
        if (e2 >-dx) { err -= dy; x0 += sx; }
        if (e2 < dy) { err += dx; y0 += sy; }
    }
}

void reshape(int width, int height) {
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glViewport(0,0,width,height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, width, 0, height);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void renderScene(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPointSize(1.0);
    glBegin(GL_POINTS);
    for(int i = 0; i < 500; i++) {
        for(int j = 0; j < 500; j++) {
            if(frameBuffer[i][j][0] > 0.005) {
                glColor3f(frameBuffer[i][j][0], frameBuffer[i][j][1], frameBuffer[i][j][2]);
                glVertex2i(i, j);
            }
        }
    }
    glEnd();
    glutSwapBuffers();
}

void count() {
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    verticesCount = 0;
    facesCount = 0;

    fp = fopen("/home/andres/viewing3d/config/poly.tri", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    while (getline(&line, &len, fp) != -1) {
        if(line[0] == 'v') {
            verticesCount++;
        } else if(line[0] == 'f') {
            facesCount++;
        }
    }
    fclose(fp);
    if (line)
        free(line);
}

void readVertex(vertex *vertices) {
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    char *vertexDefinition;
    int i = 0;

    fp = fopen("/home/andres/viewing3d/config/poly.tri", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    while (getline(&line, &len, fp) != -1) {
        if(line[0] == 'v') {
            vertexDefinition = strtok(line, " "); //ignore v
            vertexDefinition = strtok(NULL, " "); //get x
            vertices[i].x = atof(vertexDefinition);
            vertexDefinition = strtok(NULL, " "); //get y
            vertices[i].y = atof(vertexDefinition);
            vertexDefinition = strtok(NULL, " "); //get z
            vertices[i].z = atof(vertexDefinition);
            i++;
        }
    }
    fclose(fp);
    if (line)
        free(line);
}

void readFaces(polygon *faces) {
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    char *faceDefinition;
    char *vertex1;
    char *vertex2;
    char *vertex3;
    int i = 0;

    fp = fopen("/home/andres/viewing3d/config/poly.tri", "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);

    while (getline(&line, &len, fp) != -1) {
        if(line[0] == 'f') {
            faceDefinition = strtok(line, " "); //ignore f
            vertex1 = strtok(NULL, " "); //get 1
            vertex2 = strtok(NULL, " "); //get 2
            vertex3 = strtok(NULL, " "); //get 3
            if(strstr(vertex1, "/")){
                faces[i].vertices[0] = atoi(strtok(vertex1, "/"));
            } else {
                faces[i].vertices[0] = atoi(vertex1);
            }
            if(strstr(vertex1, "/")){
                faces[i].vertices[1] = atoi(strtok(vertex2, "/"));
            } else {
                faces[i].vertices[1] = atoi(vertex2);
            }
            if(strstr(vertex1, "/")){
                faces[i].vertices[2] = atoi(strtok(vertex3, "/"));
            } else {
                faces[i].vertices[2] = atoi(vertex3);
            }
            if(vertex1 == NULL || vertex2 == NULL || vertex3 == NULL) {
            }
            i++;
        }
    }
    fclose(fp);
    if (line)
        free(line);
}

int clipt(long double denom, long double num, long double *te, long double *tl){
    long double t;
    if(denom > 0) {
        t = num/denom;
        if(t > *tl) {
            return 0;
        } else if(t > *te) {
            *te = t;
        }
    } else if(denom < 0) {
        t = num/ denom;
        if(t < *te) {
            return 0;
        } else {
            *tl = t;
        }
    } else if(num > 0) {
        return 0;
    }
    return 1;
}

void clip3d(long double *x0, long double *y0, long double *z0,
            long double *x1, long double *y1, long double *z1,
            long double *zmin, int *accept) {
    long double tmin = 0.0;
    long double tmax = 1.0;
    long double dx = *x1 - *x0, dz = *z1 - *z0;
    *accept = 0;
    if(clipt(-dx -dz, *x0 + *z0, &tmin, &tmax)) {
        if(clipt(dx-dz, -*x0 + *z0, &tmin, &tmax)){
            long double dy = *y1 - *y0;
            if(clipt(dy-dz, -*y0 + *z0, &tmin, &tmax)){
                if(clipt(-dy-dz, *y0 + *z0, &tmin, &tmax)) {
                    if(clipt(-dz, *z0 - *zmin, &tmin, &tmax)) {
                        if(clipt(dz, -*z0 -1, &tmin, &tmax)) {
                            *accept = 1;
                            if(tmax < 1.0) {
                                *x1 = *x0 + tmax * dx;
                                *y1 = *y0 + tmax * dy;
                                *z1 = *z0 + tmax * dz;
                            }
                            if(tmin > 0.0) {
                                *x0 += tmin * dx;
                                *y0 += tmin * dy;
                                *z0 += tmin * dz;
                            }
                        }
                    }
                }
            }
        }
    }
}

void calculateVertex(scene_settings settings) {
    clearBuffer();
    long double b = settings.b;
    long double f = settings.f;

    matrix vrp = initVector(3);
    vrp.values[0][0] = settings.vrp[0];
    vrp.values[0][1] = settings.vrp[1];
    vrp.values[0][2] = settings.vrp[2];

    matrix vpn = initVector(3);
    vpn.values[0][0] = settings.vpn[0];
    vpn.values[0][1] = settings.vpn[1];
    vpn.values[0][2] = settings.vpn[2];

    matrix vup = initVector(3);
    vup.values[0][0] = settings.vup[0];
    vup.values[0][1] = settings.vup[1];
    vup.values[0][2] = settings.vup[2];
    vup = normalizeVector(vup);

    matrix prp = initVector(3);
    prp.values[0][0] = settings.prp[0];
    prp.values[0][1] = settings.prp[1];
    prp.values[0][2] = settings.prp[2];
    matrix prpH = getHomogeneousVector(prp);

    long double window[] = {settings.window[0], settings.window[1],
                            settings.window[2], settings.window[3]};
    long double viewport[] = {settings.viewport[0], settings.viewport[1],
                              settings.viewport[2], settings.viewport[3],
                              settings.viewport[4], settings.viewport[5]};

    matrix rz = normalizeVector(vpn);
    matrix rx = normalizeVector(crossProduct(vup, rz));
    matrix ry = crossProduct(rz, rx);

    matrix r = initMatrix(4 , 4);

    for(int i = 0; i < rx.width; i++) {
        r.values[0][i] = rx.values[0][i];
        r.values[1][i] = ry.values[0][i];
        r.values[2][i] = rz.values[0][i];
    }
    r.values[3][3] = 1;

    matrix cw = initVector(4);
    cw.values[0][0] = (window[1] + window[0])/2.0;
    cw.values[0][1] = (window[3] + window[2])/2.0;
    cw.values[0][2] = 0;
    cw.values[0][3] = 1;
    matrix dop = substract(cw, prpH);

    matrix sh = dopShearMatrix(dop);

    matrix weird = initMatrix(4, 1);
    weird.values[3][0] = 1;

    matrix temp = multiplyMatrixByMatrix(sh, translateOriginMatrix(prp));
    matrix vrp2 = multiplyMatrixByMatrix(temp, weird);

    matrix sper = initVector(4);

    sper.values[0][0] = 2*vrp2.values[2][0]/((window[1] - window[0]) * (vrp2.values[2][0] + b));
    sper.values[0][1] = 2*vrp2.values[2][0]/((window[3] - window[2]) * (vrp2.values[2][0] + b));
    sper.values[0][2] = -1/(vrp2.values[2][0] + b);
    sper = scaleMatrix(sper);

    matrix nper = multiplyMatrixByMatrix(sper, sh);
    nper = multiplyMatrixByMatrix(nper, translateOriginMatrix(prp));
    nper = multiplyMatrixByMatrix(nper, r);
    nper = multiplyMatrixByMatrix(nper, translateOriginMatrix(vrp));

    long double zmin = -((vrp2.values[0][2]+f)/(vrp2.values[0][2]+b));

    matrix mper = identityMatrix(4);
    mper.values[3][2] = -1;
    mper.values[3][3] = 0;

    matrix svv3dvTemp = initVector(3);
    svv3dvTemp.values[0][0] = (viewport[1] - viewport[0])/2;
    svv3dvTemp.values[0][1] = (viewport[3] - viewport[2])/2;
    svv3dvTemp.values[0][2] = (viewport[5] - viewport[4])/2;
    matrix svv3dv = scaleMatrix(svv3dvTemp);

    matrix translateViewPort = initVector(3);
    translateViewPort.values[0][0] = viewport[0];
    translateViewPort.values[0][1] = viewport[2];
    translateViewPort.values[0][2] = viewport[4];
    matrix translateViewPortMatrix = translate(translateViewPort);

    matrix translateCorner = initVector(3);
    translateCorner.values[0][0] = 1;
    translateCorner.values[0][1] = 1;
    translateCorner.values[0][2] = 1;
    matrix translateCornerMatrix = translate(translateCorner);

    matrix mvv3dv = multiplyMatrixByMatrix(translateViewPortMatrix, svv3dv);
    mvv3dv = multiplyMatrixByMatrix(mvv3dv, translateCornerMatrix);
    printMatrix(mvv3dv, "mvv3dv");
    for(int i = 0; i < verticesCount; i++) {
        matrix step_2 = applyTransformation(vertices[i], nper);
        transformedVertices[i].x = step_2.values[0][0];
        transformedVertices[i].y = step_2.values[0][1];
        transformedVertices[i].z = step_2.values[0][2];
    }
    for(int i = 0; i < facesCount; i++){
        vertex vertex1 = transformedVertices[faces[i].vertices[0]-1];
        vertex vertex2 = transformedVertices[faces[i].vertices[1]-1];
        vertex vertex3 = transformedVertices[faces[i].vertices[2]-1];

        long double *normal = getNormal(vertex1, vertex2, vertex3);

        transformedVertices[faces[i].vertices[0]-1].avg_normal[0] += normal[0];
        transformedVertices[faces[i].vertices[1]-1].avg_normal[0] += normal[0];
        transformedVertices[faces[i].vertices[2]-1].avg_normal[0] += normal[0];

        transformedVertices[faces[i].vertices[0]-1].avg_normal[1] += normal[1];
        transformedVertices[faces[i].vertices[1]-1].avg_normal[1] += normal[1];
        transformedVertices[faces[i].vertices[2]-1].avg_normal[1] += normal[1];

        transformedVertices[faces[i].vertices[0]-1].avg_normal[2] += normal[2];
        transformedVertices[faces[i].vertices[1]-1].avg_normal[2] += normal[2];
        transformedVertices[faces[i].vertices[2]-1].avg_normal[2] += normal[2];

        transformedVertices[faces[i].vertices[0]-1].cantFaces++;
        transformedVertices[faces[i].vertices[1]-1].cantFaces++;
        transformedVertices[faces[i].vertices[2]-1].cantFaces++;

        int accept1 = 0;
        int accept2 = 0;
        int accept3 = 0;
        clip3d(&vertex1.x, &vertex1.y, &vertex1.z,
               &vertex2.x, &vertex2.y, &vertex2.z,
               &zmin, &accept1);
        clip3d(&vertex2.x, &vertex2.y, &vertex2.z,
               &vertex3.x, &vertex3.y, &vertex3.z,
               &zmin, &accept2);
        clip3d(&vertex3.x, &vertex3.y, &vertex3.z,
               &vertex1.x, &vertex1.y, &vertex1.z,
               &zmin, &accept3);

    }
    for(int i = 0; i < verticesCount; i++) {
        transformedVertices[i] = applyAll(transformedVertices[i], mper, mvv3dv);

    }

    for(int i = 0; i < facesCount; i++) {
        vertex vertex1 = transformedVertices[faces[i].vertices[0]-1];
        vertex vertex2 = transformedVertices[faces[i].vertices[1]-1];
        vertex vertex3 = transformedVertices[faces[i].vertices[2]-1];

        vertex1.avg_normal[0] = vertex1.avg_normal[0]/(long double) vertex1.cantFaces;
        vertex1.avg_normal[1] = vertex1.avg_normal[1]/(long double) vertex1.cantFaces;
        vertex1.avg_normal[2] = vertex1.avg_normal[2]/(long double) vertex1.cantFaces;

        vertex2.avg_normal[0] = vertex2.avg_normal[0]/(long double) vertex2.cantFaces;
        vertex2.avg_normal[1] = vertex2.avg_normal[1]/(long double) vertex2.cantFaces;
        vertex2.avg_normal[2] = vertex2.avg_normal[2]/(long double) vertex2.cantFaces;

        vertex3.avg_normal[0] = vertex3.avg_normal[0]/(long double) vertex3.cantFaces;
        vertex3.avg_normal[1] = vertex3.avg_normal[1]/(long double) vertex3.cantFaces;
        vertex3.avg_normal[2] = vertex3.avg_normal[2]/(long double) vertex3.cantFaces;

        drawTriangle(vertex1, vertex2, vertex3);
    }
}


void keyboardInput(unsigned char key, int x, int y) {
    switch(key) {
        case '1':
            type = 1;
            break;
        case '0':
            type = 0;
            break;
        case 'q':
            settings.vpn[0] = settings.vpn[0] * cos(0.1) + settings.vpn[2] * sin(0.1);
            settings.vpn[2] = -settings.vpn[0] * sin(0.1) + settings.vpn[2] * cos(0.1);
            settings.vrp[0] = settings.vrp[0] * cos(0.1) + settings.vrp[2] * sin(0.1);
            settings.vrp[2] = -settings.vrp[0] * sin(0.1) + settings.vrp[2] * cos(0.1);
            break;
        case 'w':
            settings.vpn[0] = settings.vpn[0] * cos(-0.1) + settings.vpn[2] * sin(-0.1);
            settings.vpn[2] = -settings.vpn[0] * sin(-0.1) + settings.vpn[2] * cos(-0.1);
            settings.vrp[0] = settings.vrp[0] * cos(-0.1) + settings.vrp[2] * sin(-0.1);
            settings.vrp[2] = -settings.vrp[0] * sin(-0.1) + settings.vrp[2] * cos(-0.1);
            break;
        case 'y':
            settings.vrp[2] -= 1;
            break;
        case 'u':
            settings.vrp[2] += 1;
            break;
        case 'a':
            settings.vup[0] = settings.vup[0] * cos(0.1) - settings.vup[1] * sin(0.1);
            settings.vup[1] = settings.vup[0] * sin(0.1) + settings.vup[1] * cos(0.1);
            break;
        case 's':
            settings.vup[0] = settings.vup[0] * cos(-0.1) - settings.vup[1] * sin(-0.1);
            settings.vup[1] = settings.vup[0] * sin(-0.1) + settings.vup[1] * cos(-0.1);
            break;
        case 'z':
            settings.vpn[0] = 0;
            settings.vpn[1] = 0;
            settings.vpn[2] = 1;
            settings.vrp[0] = 0;
            settings.vrp[1] = 0;
            settings.vrp[2] = 50;
            settings.vup[0] = 0;
            settings.vup[1] = 1;
            settings.vup[2] = 0;
            break;
        default:
            break;
    }
    calculateVertex(settings);
    glutPostRedisplay();
}

int main(int argc, char **argv) {
    count();
    vertices = malloc(verticesCount * sizeof(vertex));
    transformedVertices = malloc(verticesCount * sizeof(vertex));
    faces = malloc(facesCount * sizeof(polygon));
    readVertex(vertices);
    readFaces(faces);

    settings.b = -1;
    settings.f = 1;
    long double vrp[] = {0,0,50};
    settings.vrp = createVector(vrp, 3);
    long double vpn[] = {0,0,1};
    settings.vpn = createVector(vpn, 3);
    long double vup[] = {0,1,0};
    settings.vup = createVector(vup, 3);
    long double prp[] = {0,0,50};
    settings.prp = createVector(prp, 3);

    long double window[] = {-50, 50, -50, 50};
    settings.window = createVector(window, 4);
    long double viewport[] = {0, 500, 0, 500, 0, 500};
    settings.viewport = createVector(viewport, 6);


    calculateVertex(settings);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(settings.viewport[1],settings.viewport[3]);
    glutCreateWindow("Drawing Polygons");
    glutDisplayFunc(renderScene);
    glutKeyboardFunc(keyboardInput);
    glutReshapeFunc(reshape);

    glutMainLoop();

    return 0;
}