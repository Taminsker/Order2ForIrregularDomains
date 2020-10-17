/** @file interface.h */

#ifndef INTERFACE_H
#define INTERFACE_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"
#include "../Data/datatypedefinitions.h"

#include "../Tools/matrix.h"
#include "../Tools/impose.h"
#include "../Tools/solver.h"
#include "../Tools/border.h"

#include "../Toolbox/differencefinite.h"
#include "../Toolbox/toolbox.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
  Énumérateur TVD
  */
typedef enum {
    TVD_RK2, // TVD- Runge Kutta d'ordre 2 en temps
    TVD_RK3, // TVD- Runge Kutta d'ordre 3 en temps
    TVD_RK4  // TVD- Runge Kutta d'ordre 4 en temps
} TVD;

/**
 * @brief Itérateur de phi sur le maillage avec le champs vectoriel W.
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param W un champs vectoriel sur le domaine.
 * @param dt le pas de temps à utiliser.
 * @param phi_n le vecteur des valeurs de phi au temps précédent.
 * @param tvd le type réprésentant l'ordre du schéma TVD-RK à utiliser.
 * @return Un vecteur représentant au temps t+dt
 */
Vector IteratePhi (Mesh* mesh, std::vector<Point*>* W, double dt, Vector* phi_n, TVD tvd = TVD_RK4);

/**
 * @brief Calcul l'approximation par WENO ordre 5 <W | Grad u>.
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param u un pointeur vers un vecteur représentant u dont il faut calculer le gradient.
 * @param W le champs vectoriel utilisé.
 * @return Un vecteur de valeur L(u) = -<W | Grad u>.
 */
Vector Weno(Mesh* mesh, Vector* u, std::vector<Point*>* W);

/**
 * @brief Réinitialise phi à l'aide de l'équation de Hamilton-Jacobi associé à la redistanciation.
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param phi le vecteur des valeurs de phi à réinitialiser.
 * @param itMax nombre de boucle temporelle à effectuer, par défaut 10.
 * @param tvd le type réprésentant l'ordre du schéma TVD-RK à utiliser.
 * @return Un vecteur représentant phi redistanciée.
 */
Vector ReInitPhi (Mesh* mesh, Vector* phi, int itMax = 10, TVD tvd = TVD_RK4);

/**
 * @brief Calcul de l'opérateur Hamiltonien H de l'équation de redistanciation à l'aide d'un flux de Godunov.
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param phi le vecteur sur lequel il faut appliquer le schéma.
 * @param phi_0 le vecteur utilisé pour calculer la fonction de lissage S.
 * @return Un vecteur de valeur H_G(phi) = S(phi_0) (1 - ||phi||) à l'aide d'un flux de Godunov.
 */
Vector Hamiltonien (Mesh* mesh, Vector* phi, Vector* phi_0);

/** @} */

#endif // INTERFACE_H
