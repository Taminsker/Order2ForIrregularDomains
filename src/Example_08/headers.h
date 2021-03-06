#ifndef HEADERS_H
#define HEADERS_H

#include "../O2FID/O2FID.h"

typedef enum {
    BackwardEuler1,
    BackwardEuler2,
    CrankNicolson,
} STYLE;

Point phi (double theta, double t = 0); // fonction levelset
Point phigrad (double theta, double t = 0);

double f (Point a, double t = 0.); // fonction de second membre
double u (Point a, double t = 0.);

void GetMatrix(Mesh* mesh, double dt, STYLE style, Matrix* Ar, Matrix* Al);

void Make(Point* minima, Point* extrema, double time_final, STYLE style, std::vector<int> listNx, std::vector<int> listNy, std::vector<int> listNz);

#endif // HEADERS_H
