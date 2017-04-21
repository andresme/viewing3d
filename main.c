#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <string.h>
#include "util/math/mathUtils.h"
#include "struct/vertex.h"

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

vertex *vertices;
polygon *faces;

void reshape(int width, int height) {
    glViewport(0,0,width,height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, width, 0, height);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}


void renderScene(void) {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLineWidth(0.5);
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    for(int i = 0; i < facesCount; i++) {

        vertex vertex1 = vertices[faces[i].vertices[0]-1];
        vertex vertex2 = vertices[faces[i].vertices[1]-1];
        vertex vertex3 = vertices[faces[i].vertices[2]-1];

        glVertex2f(vertex1.x, vertex1.y);
        glVertex2f(vertex2.x, vertex2.y);

        glVertex2f(vertex2.x, vertex2.y);
        glVertex2f(vertex3.x, vertex3.y);

        glVertex2f(vertex3.x, vertex3.y);
        glVertex2f(vertex1.x, vertex1.y);
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

    fp = fopen("/home/andres/CG/Polygons/config/poly.tri", "r");
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

    fp = fopen("/home/andres/CG/Polygons/config/poly.tri", "r");
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

    fp = fopen("/home/andres/CG/Polygons/config/poly.tri", "r");
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
                printf("something is wrong");
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
            printf("t: %Lf\n", t);
            *te = t;
            printf("te: %Lf\n", *te);
        }
    } else if(denom < 0) {
        t = num/ denom;
        if(t < *te) {
            return 0;
        } else {
            printf("t: %Lf\n", t);
            *tl = t;
            printf("tl: %Lf\n", *tl);
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
        printf("here1\n");
        if(clipt(dx-dz, -*x0 + *z0, &tmin, &tmax)){
            printf("here2\n");
            long double dy = *y1 - *y0;
            if(clipt(dy-dz, -*y0 + *z0, &tmin, &tmax)){
                printf("here3\n");
                if(clipt(-dy-dz, *y0 + *z0, &tmin, &tmax)) {
                    printf("here4\n");
                    if(clipt(-dz, *z0 - *zmin, &tmin, &tmax)) {
                        printf("here5\n");
                        if(clipt(dz, -*z0 -1, &tmin, &tmax)) {
                            printf("tmax: %Lf, tmin: %Lf\n", tmax, tmin);
                            *accept = 1;
                            if(tmax < 1.0) {
                                printf("x: %Lf, y: %Lf, z: %Lf\n", *x1, *y1, *z1);
                                *x1 = *x0 + tmax * dx;
                                *y1 = *y0 + tmax * dy;
                                *z1 = *z0 + tmax * dz;
                                printf("x: %Lf, y: %Lf, z: %Lf\n", *x1, *y1, *z1);
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


int main(int argc, char **argv) {
    count();
    vertices = malloc(verticesCount * sizeof(vertex));
    faces = malloc(facesCount * sizeof(polygon));
    readVertex(vertices);
    readFaces(faces);

    long double b = -23;
    long double f = 24;
    matrix vrp = initVector(3);
    vrp.values[0][0] = 0;
    vrp.values[0][1] = 0;
    vrp.values[0][2] = 54;
    matrix vrpH = getHomogeneousVector(vrp);

    matrix vpn = initVector(3);
    vpn.values[0][0] = 0;
    vpn.values[0][1] = 0;
    vpn.values[0][2] = 1;
    matrix vpnH = getHomogeneousVector(vpn);

    matrix vup = initVector(3);
    vup.values[0][0] = 0;
    vup.values[0][1] = 1;
    vup.values[0][2] = 0;
    matrix vupH = getHomogeneousVector(vup);

    matrix prp = initVector(3);
    prp.values[0][0] = 8;
    prp.values[0][1] = 6;
    prp.values[0][2] = 30;
    matrix prpH = getHomogeneousVector(prp);

    long double window[] = {-1, 1, -1, 1};
    long double viewport[] = {0, 500, 0, 500, 0, 500};

    printMatrix(vrpH, "vrpH");
    printMatrix(translateOriginMatrix(vrp), "translate matrix");
    matrix vrpOrigin = transpose(multiplyMatrixByMatrix(translateOriginMatrix(vrp), transpose(vrpH)));
    printMatrix(vrpOrigin,"origin vrp");

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

    printMatrix(rx, "rx");
    printMatrix(ry, "ry");
    printMatrix(rz, "rz");
    printMatrix(r, "r");

    matrix prpOrigin = transpose(multiplyMatrixByMatrix(translateOriginMatrix(prp), transpose(prpH)));
    printMatrix(prpOrigin, "prp Origin");

    matrix cw = initVector(4);
    cw.values[0][0] = (window[1] + window[0])/2.0;
    cw.values[0][1] = (window[3] + window[2])/2.0;
    cw.values[0][2] = 0;
    cw.values[0][3] = 1;
    printMatrix(cw, "cw");
    matrix dop = substract(cw, prpH);
    printMatrix(dop, "dop");

    matrix sh = dopShearMatrix(dop);
    printMatrix(sh, "sh");

    matrix weird = initMatrix(4, 1);
    weird.values[3][0] = 1;

    matrix temp = multiplyMatrixByMatrix(sh, translateOriginMatrix(prp));
    matrix vrp2 = multiplyMatrixByMatrix(temp, weird);
    printMatrix(vrp2, "vrp'");

    matrix sper = initVector(4);

    sper.values[0][0] = 2*vrp2.values[2][0]/((window[1] - window[0]) * (vrp2.values[2][0] + b));
    sper.values[0][1] = 2*vrp2.values[2][0]/((window[3] - window[2]) * (vrp2.values[2][0] + b));
    sper.values[0][2] = -1/(vrp2.values[2][0] + b);
    sper = scaleMatrix(sper);
    printMatrix(sper, "Sper");

    matrix nper = multiplyMatrixByMatrix(sper, sh);
    printMatrix(nper, "nper1");
    nper = multiplyMatrixByMatrix(nper, translateOriginMatrix(prp));
    printMatrix(nper, "nper2");
    nper = multiplyMatrixByMatrix(nper, r);
    printMatrix(nper, "nper3");
    nper = multiplyMatrixByMatrix(nper, translateOriginMatrix(vrp));
    printMatrix(nper, "nper4");

    long double zproj = -(vrp2.values[0][2]/(vrp2.values[0][2]+b));
    long double zmin = -((vrp2.values[0][2]+f)/(vrp2.values[0][2]+b));
    long double zmax = -((vrp2.values[0][2]+b)/(vrp2.values[0][2]+b));

    printf("zproj: %Lf, zmin: %Lf, zmax: %Lf\n", zproj, zmin, zmax);

    matrix mper = identityMatrix(4);
    mper.values[3][2] = -1;
    mper.values[3][3] = 0;
    printMatrix(mper, "mper");
    

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
    printMatrix(mvv3dv, "mvv3dv1");
    mvv3dv = multiplyMatrixByMatrix(mvv3dv, translateCornerMatrix);
    printMatrix(mvv3dv, "mvv3dv2");

    for(int i = 0; i < verticesCount; i++) {
        matrix step_2 = applyTransformation(vertices[i], nper);
        vertices[i].x = step_2.values[0][0];
        vertices[i].y = step_2.values[0][1];
        vertices[i].z = step_2.values[0][2];
    }
    for(int i = 0; i < facesCount; i++){
        vertex vertex1 = vertices[faces[i].vertices[0]-1];
        vertex vertex2 = vertices[faces[i].vertices[1]-1];
        vertex vertex3 = vertices[faces[i].vertices[2]-1];
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
        vertices[i] = applyAll(vertices[i], mper, mvv3dv);
    }

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(viewport[1],viewport[3]);
    glutCreateWindow("Drawing Polygons");
    glutDisplayFunc(renderScene);
    glutReshapeFunc(reshape);

    glutMainLoop();

    return 0;
}


