/** \file border.h */
#ifndef BORDER_H
#define BORDER_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
 * @brief Cette fonction trouve les noeuds qui sont supposés définir la frontière \f$\Omega\f$ et elle les rajoute au maillage.
 * @param mesh un pointeur vers un objet Mesh
 * @param f un pointeur vers une fonction f
 * @return Vecteur les indices des points sur cette frontière
 */
std::vector <int> MakeListOfIndexPoints (Mesh * mesh, double (*f) (Point, double));

/*! @} End of Doxygen Groups*/

#endif // BORDER_H
