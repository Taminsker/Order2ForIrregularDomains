/** \file border.h */
#ifndef BORDER_H
#define BORDER_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
 * @brief Cette fonction trouve les noeuds qui sont supposés définir la frontière \f$\Omega\f$ et elle les rajoute au maillage.
 * @param mesh un pointeur vers un objet Mesh
 * @param phi_list un pointeur vers un vecteur des valeurs de phi
 * @return Vecteur les indices des points sur cette frontière
 */
std::vector <int> MakeBorderPoints (Mesh * mesh, Vector * phi_list);



/*! @} End of Doxygen Groups*/

#endif // BORDER_H
