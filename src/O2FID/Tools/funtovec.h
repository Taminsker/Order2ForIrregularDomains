/** \file funtovec.h */
#ifndef FUNTOVEC_H
#define FUNTOVEC_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"

#include <vector>
#include <math.h>

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
 * @brief Cette fonction prend en argument un pointeur vers un object Mesh,
          une fonction qui retourne un double et qui prend en arguments un point
          de la classe Point et un double représentant le temps .
 * @param mesh pointeur vers un maillage.
 * @param f pointeur vers une fonction définit sur des coordonnées spatiales et une cooordonnée temporelle t.
 * @param t la coordonnée temporelle utilisée pour construire le vecteur.
 * @return Vector vecteur sur le maillage.
 * @see La classe Point et la classe Mesh.
 */
Vector FunToVec (Mesh * mesh, double (*f) (Point, double), double t = 0.);

/**
 * @brief Cette fonction prend en argument un pointeur vers un object Mesh et un simple double.
 * @param mesh pointeur vers un maillage.
 * @param value valeur qui initialise tout le vecteur à cette valeur.
 * @return Vector vecteur sur le maillage initialisé à la valeur value.
 * @see La classe Point et la classe Mesh.
 *
 */
Vector FunToVec (Mesh * mesh, double value = 0.);


typedef struct {
    Vector* out;
    std::vector<Point*>* position;
    std::vector<Point*>* velocity;
    std::vector <bool>* reached;
    size_t size;
    double ratio;
} WagonF2V;

/**
 * @brief Cette fonction prend en argument un pointeur vers un object Mesh,
          une fonction qui retourne un double et qui prend en arguments un point
          de la classe Point et un double représentant le temps .
 * @param mesh pointeur vers un maillage.
 * @param f pointeur vers une fonction définit sur une coordonnée angulaire theta et une cooordonnée temporelle t.
 * @param t la coordonnée temporelle utilisée pour construire le vecteur.
 * @return Vector vecteur sur le maillage.
 * @see La classe Point et la classe Mesh.
 */
Vector FunToVec (Mesh * mesh,
                 Point (*f) (double, double),
                 Point (*fprim) (double, double),
                 size_t numOfPtsBorder = 2000,
                 double t = 0,
                 double a = 0,
                 double b = 2 * M_PI);


std::vector<std::vector <int>> Build (int i, std::pair<int, int> j, int k);

std::vector<std::vector <int>> Build (std::pair<int, int> i, int j, int k);


int Remainder (int dividend, int divisor);

/*! @} End of Doxygen Groups*/

#endif // FUNTOVEC_H
