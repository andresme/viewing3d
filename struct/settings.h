#ifndef POLYGONS_SETTINGS_H
#define POLYGONS_SETTINGS_H


typedef struct scene_settings {
    long double b;
    long double f;
    long double *vrp;
    long double *vpn;
    long double *prp;
    long double *vup;
    long double *window;
    long double *viewport;
}scene_settings;

long double* createVector(long double vector[], int size);

#endif //POLYGONS_SETTINGS_H
