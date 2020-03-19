/** @file point.h */

#ifndef POINT_H
#define POINT_H

#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <initializer_list>

#define SPACE std::left << std::setw(12)

/**
 * @brief Rappel de l'existence de la classe Cell (le header n'est pas inclu directement car il y avait un risque d'inclusion réciproque au moment où je codais cette classe).
 */
class Cell;

/**
 * @brief Enumerateur pour la localisation des points.
 */
typedef enum
{
    ON_BORDER_OMEGA,
    ON_DOMAIN_INTERN_OMEGA,
    ON_DOMAIN_EXTERN_OMEGA
} POINT_LOCATION;

/**
 * @brief La classe Point implémente les informations nécessaires pour définir un point.
 * Cette classe permet des points 3D via plusieurs constructeurs. L'opérateur d'égalité est aussi surchargé ainsi que l'opérateur d'affichage.
 */
class Point
{

public:
    /**
     * @brief Coordonnée x du point.
     */
    double x = 0.;

    /**
     * @brief Coordonnée z du point.
     */
    double y = 0.;

    /**
     * @brief Coordonnée y du point.
     */
    double z = 0.;

    /**
     * @brief Créer un point (0, 0, 0).
     */
    Point ();

    /**
     * @brief Créer une copie du point entré p.
     * @param p Point
     */
    Point (const Point &p);

    /**
     * @brief Créer un point (a, b, c).
     * @param a double.
     * @param b double.
     * @param c double.
     */
    Point (double a, double b = 0., double c = 0.);

    /**
     * @brief Destructeur
     */
    ~Point ();

    /**
     * @brief Définition de l'indice global du point
     * @param index un entier int
     * @return rien
     */
    void SetGlobalIndex (int index);

    /**
     * @brief Retourne l'indice global du point dans le maillage auquel il appartient
     * @return index l'indice global
     */
    int GetGlobalIndex ();

    /**
     * @brief Définition de la localisation du point. Voir l'énumération POINT_LOCATION
     * @param loc tag de type POINT_LOCATION
     * @return rien
     */
    void SetLocate (POINT_LOCATION loc);

    /**
     * @brief Retourne le tag du point.
     * @return tag de type POINT_LOCATION
     */
    POINT_LOCATION GetLocate () const;

    /**
     * @brief Nettoie la liste des points voisins de ce point
     * @return rien
     */
    void ClearListNeighboors ();

    /**
     * @brief Ajoute de le point p à la liste de voisins de ce point
     * @param p point de type Point
     * @return rien
     */
    void AddPointNeighbour (Point p);

    /**
     * @brief Retourne un vecteur de points voisins de ce point
     * @return vecteur de point
     */
    std::vector <Point> GetListNeighbours ();

    /**
     * @brief operator == : operateur de test d'égalité sur les 3 coordonnées.
     * @param p Point
     */
    bool operator== (const Point &p);

    /**
     * @brief Addtition coordonnée par coordonnée
     * @param a un point de type Point
     * @return Point
     */
    Point operator+ (const Point & a);

    /**
     * @brief Soustraction coordonnée par coordonnée
     * @param a un point de type Point
     * @return Point
     */
    Point operator- (const Point & a);

    /**
     * @brief Multiplication par un double de toutes les coordonnées
     * @param v un double
     * @return Point
     */
    Point operator* (double v);

    /**
     * @brief Division par un double de toutes les coordonnées (si différent de zéro)
     * @param v un double
     * @return Point
     */
    Point operator/ (double v);

    /**
     * @brief Multiplication par un int de toutes les coordonnées
     * @param v un double
     * @return Point
     */
    Point operator* (int v);

    /**
     * @brief operator << : utiliser pour std::cout << p << std::endl, similaire à printf ("%f, %f, %f", p.x, p.y, p.z).
     * @param p Point
     * @param out flux
     * @return flux
     */
    friend std::ostream & operator<< (std::ostream &out, const Point &p);

protected:
    /**
     * @brief Liste des cellules connectées en ce point
     */
    std::vector <Cell *> m_cells;

    /**
     * @brief Tag de localisation du point
     */
    POINT_LOCATION m_locate;

    /**
     * @brief Indice global du point dans le maillage
     */
    int m_globalIndex;

    /**
    * @brief Liste des points voisins
    */
    std::vector<Point> m_listNeighbours;
};

/**
 * @brief Opérateur d'affichage (surcharge des doubles chevrons)
 * operator << : utiliser pour std::cout << p << std::endl, similaire à printf ("%f, %f, %f", p.x, p.y, p.z).
 * @param p Point
 * @param out flux
 * @return flux
 */
std::ostream & operator<< (std::ostream &out, const Point &p);

#endif // POINT_H
