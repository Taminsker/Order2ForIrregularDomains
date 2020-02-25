/** \file funtovec.h */
#ifndef FUNTOVEC_H
#define FUNTOVEC_H

#include "../Data/data.h"
#include "../Mesh/mesh.h"

/**
 * @brief La classe FunToVec sert à créer un vecteur sur un maillage
 * Cette classe sert notamment à créer des vecteurs sur le maillage donné à la création de la classe.
 * Plusieurs fonctions 'Make' sont disponibles :
 * - fonction prenant en argument un simple point de la classe Point
 * - fonction prenant en arguments un point et un double représentant le temps
 * - ou simplement une fonction retournant un vecteur initialisé à une certaine valeur passée en argument
 * @see La classe Point
 */
class FunToVec
{
public:
    /**
     * @brief Constructeur prenant en argument un pointeur vers un maillage.
     * @arg maillage pointeur vers un maillage de type Mesh.
     * @see La classe Mesh.
     */
    FunToVec (Mesh * mesh);

    /**
     * @brief Destructeur de cette classe
     */
    ~FunToVec ();

    /**
     * @brief Cette méthode 'Make' prend en argument une fonction qui retourne un double et qui prend en argument un point de la classe Point.
     * @param f pointeur vers une fonction uniquement définit sur des coordonnées spatiales.
     * @return Vector : vecteur sur le maillage.
     * @see La classe Point.
     */
    Vector Make (double (*f) (Point));

    /**
     * @brief Cette méthode 'Make' prend en argument un double représentant le temps et une fonction qui retourne un double et qui prend en arguments un point de la classe Point et un double représentant le temps .
     * @param f pointeur vers une fonction définit sur des coordonnées spatiales et une cooordonnée temporelle t.
     * @param t la coordonnée temporelle utilisée pour construire le vecteur.
     * @return Vector vecteur sur le maillage.
     * @see La classe Point.
     */
    Vector Make (double (*f) (Point, double), double t);

    /**
     * @brief Cette méthode 'Make' prend un argument un simple double.
     * @param value valeur qui initialise tout le vecteur à cette valeur.
     * @return Vector vecteur sur le maillage.
     * @see La classe Point.
     */
    Vector Make (double value);

protected:
    /**
     * @brief Pointeur vers un maillage.
     * @see La classe Mesh.
     */
    Mesh * m_mesh;
};

#endif // FUNTOVEC_H
