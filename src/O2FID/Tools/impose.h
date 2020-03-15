/** \file impose.h */
#ifndef IMPOSE_H
#define IMPOSE_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"

#include "Eigen/SparseCore"

/*!
 *  \addtogroup Outils
 *  @{
 */


/**
 * @brief Cette fonction modifie la sparsematrix pour imposer la conditions de Dirichlet et elle retourne le vecteur des conditions imposées (voir calculs).
 * @param mesh un pointeur vers un objet Mesh
 * @param sparsematrix un pointeur vers une sparsematrix d'Eigen
 * @param g un pointeur vers une fonction g dépendant des coordonnées spatiales et temporelle
 * @param t un réel (optionnel, par défaut il vaut 0)
 * @param listIndex un vecteur d'entiers regroupant les indices des points comme définissant la frontière de \f$\Omega\f$
 * @return le vecteur des conditions imposées
 */
Vector ImposeDirichlet (Mesh &mesh,
                        Matrix * sparsematrix,
                        double (*g) (Point, double),
                        std::vector <int> listIndex,
                        double t = 0.);

/**
 * @brief Cette fonction modifie la sparsematrix pour imposer la conditions de Neumann et elle retourne le vecteur des conditions imposées (voir calculs).
 * @param mesh un pointeur vers un objet Mesh
 * @param sparsematrix un pointeur vers une sparsematrix d'Eigen
 * @param g un pointeur vers une fonction g dépendant des coordonnées spatiales et temporelle
 * @param phi un pointeur vers une fonction phi de levelset
 * @param t un réel (optionnel, par défaut il vaut 0)
 * @param listIndex un vecteur d'entiers regroupant les indices des points comme définissant la frontière de \f$\Omega\f$
 * @return le vecteur des conditions imposées
 */
Vector ImposeNeumann (Mesh &mesh,
                      Matrix * sparsematrix,
                      double (*g) (Point, double),
                      double (*phi) (Point, double),
                      std::vector <int> listIndex,
                      double t = 0.);


/*! @} End of Doxygen Groups*/

#endif // IMPOSE_H
