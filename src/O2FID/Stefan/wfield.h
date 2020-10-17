/** @file wfield.h */

#ifndef WFIELD_H
#define WFIELD_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"
#include "../Data/datatypedefinitions.h"

#include "../Tools/matrix.h"
#include "../Tools/impose.h"
#include "../Tools/solver.h"
#include "../Tools/border.h"

#include "../Toolbox/differencefinite.h"

#include "field.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
 * @brief Construit le champs de vecteur W sur le domaine.
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param phi un pointeur vers un vecteur de valeur représentant phi sur le maillage.
 * @param sol un pointeur vers un vecteur de solution numérique.
 * @param idxsBorder liste des indices de points de bords.
 * @param h0 réel associé au problème.
 * @return Un pointeur vers un objet Field.
 */
Field* GetWField (Mesh* mesh, Vector* phi, Vector* sol, std::vector<int>* idxsBorder, double h0 = 1);

/**
 * @brief Contruit le champs de vecteur W sur le bord
 * @param field pointeur vers l'objet regroupant les différents vecteurs.
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param phi un pointeur vers un vecteur de valeur représentant phi sur le maillage.
 * @param idxsBorder liste des indices de points de bords.
 * @param h0 réel associé au problème.
 * @return Un vecteur de pointeur de Point, ici les W.
 */
std::vector<Point*> BuildWOnBorder (Field* field, Mesh* mesh, Vector* phi, std::vector<int>* idxsBorder, double h0);

/**
 * @brief Étend le champs de vecteur W à tout le domaine à partir des valeurs sur le bord (résoud Lap W = 0 avec conditions de Dirichlet).
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param idxsBorder liste des indices de points de bords.
 * @param W le champs de vecteurs associé au domaine.
 */
void ExtendWToAllDomain (Mesh* mesh, std::vector<int>* idxsBorder, std::vector<Point*>* W);

/**
 * @brief Calcul de gradient associé associé au vecteur vec sur le domaine.
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param vec un pointeur vers un vecteur de valeur sur le maillage.
 * @param normalized booléen de normalisation.
 * @param order ordre associé à l'approximation @see DFStruct.
 * @return vecteur de gradient de vec en chaque point du maillage.
 */
std::vector<Point*> ComputeGradient (Mesh* mesh, Vector* vec, bool normalized, ORDERS order = ORDER_2_CENTRAL);

/**
 * @brief Calcul de la norme d'un champs de vecteurs sur le maillage.
 * @param vec le champs de vecteurs dont la norme va être calculée.
 * @return Un vecteur de valeurs réelles.
 */
Vector ComputeNorm(std::vector<Point*>* vec);

/**
 * @brief Calcul du dt optimal en fonction d'un champs vectoriel.
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param field le champs vectoriel utilisé (voir la condition de CFL).
 * @return dt optimal.
 */
double Compute_dt(Mesh* mesh, Field* field);

/** @} */

#endif // WFIELD_H
