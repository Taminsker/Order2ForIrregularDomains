/** \file impose.h */\


#ifndef IMPOSE_H
#define IMPOSE_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"
#include "../Data/datatypedefinitions.h"

#include "Eigen/SparseCore"

/*!
 *  \addtogroup Outils
 *  @{
 */


/**
 * @brief Cette fonction modifie la sparsematrix pour imposer la conditions de Dirichlet et elle retourne le vecteur des conditions imposées (voir calculs).
 * @param mesh un pointeur vers un objet Mesh
 * @param A un pointeur vers une sparsematrix d'Eigen
 * @param secondMember vecteur de second membre
 * @param g un pointeur vers une fonction g dépendant des coordonnées spatiales et temporelle, ou bien g_list le vecteur de g sur le maillage
 * @param t un réel (optionnel, par défaut il vaut 0)
 * @param listIndex un vecteur d'entiers regroupant les indices des points comme définissant la frontière de \f$\Omega\f$
 */
void ImposeDirichlet (Mesh* mesh,
                      Matrix* A,
                      Vector* secondMember,
                      double (*g) (Point, double),
                      std::vector <int> listIndex,
                      double t = 0.);



void ImposeDirichlet (Mesh* mesh,
                      Matrix* A,
                      Vector* secondMember,
                      Vector* g_list,
                      std::vector <int> listIndex,
                      double t = 0.);


/**
  * @brief Énumération sur le degré de l'interpolation utilisée pour imposer une condition de Neumann
  */
typedef enum {
    /** Interpolation de degré 1 */
    DEGRE_1 = 1,
    /** Interpolation de degré 2 */
    DEGRE_2 = 2,
    /** Interpolation de degré 3 */
    DEGRE_3 = 3,
    /** Interpolation de degré 4 */
    DEGRE_4 = 4
} INTERPOLATION_NORMAL;



/**
 * @brief Cette fonction modifie la sparsematrix pour imposer la conditions de Neumann et elle retourne le vecteur des conditions imposées (voir calculs).
 * @param mesh un pointeur vers un objet Mesh
 * @param A pointeur vers une sparsematrix d'Eigen
 * @param scMember vecteur Eigen de second membre
 * @param g un pointeur vers une fonction g dépendant des coordonnées spatiales et temporelle
 * @param phigrad un pointeur vers le gradient d'une fonction phi de levelset
 * @param listOfIndexes un vecteur d'entiers regroupant les indices des points comme définissant la frontière de \f$\Omega\f$
 * @param interpolationType type d'interpolation sélectionné
 * @param t un réel (optionnel, par défaut il vaut 0)
 */
void ImposeNeumann (  Mesh *mesh,
                      Matrix * A,
                      Vector* scMember,
                      double (*g) (Point, double),
                      Point (*phigrad) (Point, double),
                      std::vector <int> listOfIndexes,
                      INTERPOLATION_NORMAL interpolationType = DEGRE_1,
                      double t = 0.);

/**
 * @brief Cette fonction modifie la sparsematrix pour imposer la conditions de Neumann et elle retourne le vecteur des conditions imposées (voir calculs).
 * @param mesh un pointeur vers un objet Mesh
 * @param A pointeur vers une sparsematrix d'Eigen
 * @param scMember vecteur Eigen de second membre
 * @param g un pointeur vers une fonction g dépendant des coordonnées spatiales et temporelle
 * @param phi un pointeur vers une fonction phi de levelset
 * @param phigrad un pointeur vers une le gradient de la fonction phi
 * @param listOfIndexes un vecteur d'entiers regroupant les indices des points comme définissant la frontière de \f$\Omega\f$
 * @param interpolationType type d'interpolation sélectionné
 * @param t un réel (optionnel, par défaut il vaut 0)
 */
void ImposeNeumann (  Mesh *mesh,
                      Matrix * A,
                      Vector* scMember,
                      double (*g) (Point, double),
                      Point (*phi) (double, double),
                      Point (*phigrad) (double, double),
                      std::vector <int> listOfIndexes,
                      INTERPOLATION_NORMAL interpolationType = DEGRE_1,
                      double t = 0.);

void ImposeZeroDirichletExtBorder (Mesh* mesh, Matrix* A, Vector* a1, Vector* a2 = nullptr, Vector* a3 = nullptr);


/*! @} End of Doxygen Groups*/

#endif // IMPOSE_H
