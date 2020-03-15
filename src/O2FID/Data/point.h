/** @file point.h */

#ifndef POINT_H
#define POINT_H

#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>

#define SPACE std::left << std::setw(12)

class Cell;

/**
 * @brief Enumerateur pour la localisation des points.
 */
typedef enum
{
    ON_BORDER_OMEGA,
    ON_DOMAIN_OMEGA
} POINT_LOCATION;

/**
 * @brief La classe Point implémente les informations nécessaires pour définir un point.
 * Cette classe met permet des points 3D via plusieurs constructeurs. L'opérateur d'égalité est aussi surchargé ainsi que l'opérateur d'affichage.
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
     * @brief SetLocate
     * @param loc
     */
    void SetLocate (int loc);

    /**
     * @brief GetLocate
     * @return
     */
    int GetLocate () const;

    /**
     * @brief operator == : operateur de test d'égalité sur les 3 coordonnées.
     * @param p Point
     */
    bool operator== (const Point &p);

    /**
     * @brief operator +
     * @param a
     * @return
     */
    Point operator+ (const Point & a);

    /**
     * @brief operator -
     * @param a
     * @return
     */
    Point operator- (const Point & a);

    /**
     * @brief operator *
     * @param v
     * @return
     */
    Point operator* (double v);

    /**
     * @brief operator /
     * @param v
     * @return
     */
    Point operator/ (double v);

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
     * @brief m_locate
     */
    int m_locate;
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
