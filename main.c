#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <string.h>
#include "util/math/mathUtils.h"
#include "struct/vertex.h"
#include "struct/settings.h"

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
vertex *transformedVertices;
polygon *faces;

scene_settings settings;

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

        vertex vertex1 = transformedVertices[faces[i].vertices[0]-1];
        vertex vertex2 = transformedVertices[faces[i].vertices[1]-1];
        vertex vertex3 = transformedVertices[faces[i].vertices[2]-1];

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

    matrix prp = initVector(3);
    prp.values[0][0] = settings.prp[0];
    prp.values[0][1] = settings.prp[1];
    prp.values[0][2] = settings.prp[2];
    matrix prpH = getHomogeneousVector(prp);

    long double window[] = {settings.window[0], settings.window[1], settings.window[2], settings.window[3]};
    long double viewport[] = {settings.viewport[0], settings.viewport[1], settings.viewport[2], settings.viewport[3], settings.viewport[4], settings.viewport[5]};

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
}


void keyboardInput(unsigned char key, int x, int y) {
    switch(key) {
        case 'q':
            settings.vpn[0] += 0.1;
            break;
        case 'w':
            settings.vpn[0] -= 0.1;
            break;
        case 'e':
            settings.vpn[1] += 0.1;
            break;
        case 'r':
            settings.vpn[1] -= 0.1;
            break;
        case 't':
            settings.vpn[2] -= -0.1;
            break;
        case 'y':
            settings.vpn[2] -= 0.1;
            break;
        case 'a':
            settings.vrp[0] += 0.1;
            break;
        case 's':
            settings.vrp[0] -= 0.1;
            break;
        case 'd':
            settings.vrp[1] += 0.1;
            break;
        case 'f':
            settings.vrp[1] -= 0.1;
            break;
        case 'g':
            settings.vrp[2] -= -0.1;
            break;
        case 'h':
            settings.vrp[2] -= 0.1;
            break;
        case 'x':
            settings.prp[0] += 0.1;
            break;
        case 'c':
            settings.prp[0] -= 0.1;
            break;
        case 'v':
            settings.prp[1] += 0.1;
            break;
        case 'b':
            settings.prp[1] -= 0.1;
            break;
        case 'n':
            settings.prp[2] -= -0.1;
            break;
        case 'm':
            settings.prp[2] -= 0.1;
            break;

        case 'z':
            settings.vpn[0] = 0;
            settings.vpn[1] = 0;
            settings.vpn[2] = 1;
            settings.vrp[0] = 0;
            settings.vrp[1] = 0;
            settings.vrp[2] = 1;
            settings.vrp[0] = 0;
            settings.vrp[1] = 0;
            settings.vrp[2] = 50;
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
    long double vrp[] = {0,0,0};
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


