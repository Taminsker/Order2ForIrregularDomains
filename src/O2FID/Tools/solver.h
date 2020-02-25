/** \file solver.h */
#ifndef SOLVER_H
#define SOLVER_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"

/*!
 *  \addtogroup Outils
 *  @{
 */

typedef enum {
    EXPLICIT,
    IMPLICIT
} TYPE;

/**
 * @brief Cette fonction résout au choix le problème implicite avec le gradient conjugué d'Eigen ou le simple produit matrice-vecteur dans le cas explicite et elle retourne le vecteur résultant.
 * @param A une sparsematrix
 * @param b un vecteur symbolisant à la fois le second membre dans le cas implicite et la valeur au temps précédent dans le cas explicite
 * @param type un entier (valant EXPLICIT ou IMPLICIT) pour préciser le type de résolution
 * @return le vecteur résultant
 */
Vector Solve (const Matrix &A, const Vector &b, int type = EXPLICIT);

/*! @} End of Doxygen Groups*/

#endif // SOLVER_H
