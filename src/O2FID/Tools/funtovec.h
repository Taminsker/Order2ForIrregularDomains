/** \file funtovec.h */
#ifndef FUNTOVEC_H
#define FUNTOVEC_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"

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

/*! @} End of Doxygen Groups*/

#endif // FUNTOVEC_H
