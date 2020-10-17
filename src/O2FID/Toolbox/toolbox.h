/** @file toolbox.h */

#ifndef TOOLBOX_H
#define TOOLBOX_H

#define INDENT "-->\t"

#include <iostream>
#include <vector>
#include <cmath>

#include "../Data/datatypedefinitions.h"

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
 * @brief Calcul du reste de la division euclidienne
 * @param dividend le dividende
 * @param divisor le diviseur
 * @return le reste de la division
 */
int Remainder (int dividend, int divisor);

/**
 * @brief Calcul du quotient de la division euclidienne
 * @param dividend le dividende
 * @param divisor le diviseur
 * @return le quotient de la division
 */
int Quotient (int dividend, int divisor);

/**
 * @brief Calcul des ordres de convergences par rapport à un vecteur d'erreurs et de pas h.
 * @param err le vecteur des erreurs.
 * @param h le vecteur des pas d'espaces.
 * @return Un vecteur d'ordre de convergence.
 */
std::vector<double> Order (std::vector<double> err, std::vector<double> h);

/**
 * @brief Fonction Template de flux d'affichage d'un std vector
 */
template<typename T>
std::ostream & operator<< (std::ostream &out, const std::vector<T> vec)
{
    for (size_t i = 0; i < vec.size (); ++i)
        out << vec.at (i) << " " << std::flush;
    return out;
}

/**
 * @brief Fonction Template de nettoyage d'un std vector de pointeurs.
 */
class Point;
template <typename T>
std::vector<T*>* AutoClearVector(std::vector<T*>* vec)
{
    if (vec == nullptr)
        return vec;

   for (typename std::vector<T*>::iterator it =vec->begin (); it != vec->end ();)
    {
        delete *it;
        it = vec->erase(it);
    }

    return vec;
}
/**
 * @brief Fonction Template de multiplication entre un scalaire et un std vector.
 */
template <typename T>
std::vector<T> operator* (T value, std::vector<T> vec)
{
    std::vector<T> R(vec.size (), T(0));

    for (size_t i = 0; i < vec.size (); ++i)
    {
        R.at (i) = value * vec.at (i);
    }

    return R;
}

class Mesh;

/**
 * @brief Extrapole le vecteur vec sur le maillage (extrapole aux points de bords notamment).
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param vec le vecteur à étendre.
 */
void Extrapole (Mesh* mesh, Vector* vec);

/**
 * @brief Extrapole le vecteur vec sur le maillage (extrapole aux points de bords notamment).
 * @param mesh un pointeur vers un objet de type Mesh.
 * @param vec le vecteur à étendre.
 */
void Extrapole (Mesh* mesh, std::vector<Point*>* vec);

/**
 * @brief Ajoute au vecteur vec la liste des voisins de points de bords (voir fonction pour le degré de recherche), sans répétition.
 * @param mesh un pointeur vers un obket de type Mesh.
 * @param vec le vecteur d'indices de Point.
 */
void ExtendToNeighbours(Mesh* mesh, std::vector<int>* vec);

/** @} */


#endif // TOOLBOX_H
