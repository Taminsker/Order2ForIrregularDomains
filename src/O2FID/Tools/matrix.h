/** \file matrix.h */
#ifndef MATRIX_H
#define MATRIX_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
 * @brief Cette fonction retourne une sparsematrix du problème de Laplace sans préoccupation des conditions aux limites
 * @param mesh un pointeur vers un objet Mesh
 * @return Matrix du problème de Laplace
 */
Matrix BuildMatrixLaplaceEquation (Mesh * mesh);

/**
 * @brief Cette fonction retourne une sparsematrix du problème de Laplace sans préoccupation des conditions aux limites
 * @param mesh un pointeur vers un objet Mesh
 * @param dt un réel
 * @param coeff le coefficient de l'équation de la chaleur
 * @return Matrix du problème de Laplace
 */
Matrix BuildMatrixHeatEquation (Mesh * mesh, double dt, double coeff);

/*! @} End of Doxygen Groups*/

#endif // MATRIX_H
