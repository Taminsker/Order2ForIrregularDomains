/** @file point.h */

#ifndef POINT_H
#define POINT_H

#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <initializer_list>
#include <utility>

#define SPACE std::left << std::setw(12)

/**
 * @brief Rappel de l'existence de la classe Cell (le header n'est pas inclu directement car il y avait un risque d'inclusion réciproque au moment où je codais cette classe).
 */
class Cell;

/*!
 *  \addtogroup Outils
 *  @{
 */
/**
 * @brief Enumerateur pour la localisation des points.
 */
typedef enum
{
    /** Dans le domaine interne phi négative */
    ON_BORDER_OMEGA = 1,
    /** Sur le bord du domaine, phi égale à zéro */
    ON_DOMAIN_INTERN_OMEGA = 0,
    /** Dans le domaine externe phi positive */
    ON_DOMAIN_EXTERN_OMEGA = 2
} POINT_LOCATION;

/**
 * @brief Enumérateur pour le label des axes
 */
typedef enum
{
    /** Pas d'axe */
    NO_AXIS = -1,
    /** Axe des x */
    AXIS_X = 0,
    /** Axe des y */
    AXIS_Y = 1,
    /** Axe des z */
    AXIS_Z = 2
} AXIS_LABEL;

typedef enum
{
    /** Tag à gauche */
    LEFT = 0,
    /** Tag à droite */
    RIGHT = 1,
    /** Tag au dessus */
    UP = 2,
    /** Tag au dessous */
    DOWN = 3,
    /** Tag de devant */
    FRONT = 4,
    /** Tag tag de derrière */
    BACK = 5
} POSITION;

/** @}*/

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
     * @brief Affectation par égalité
     * @param p Point
     * @return
     */
    Point& operator= (const Point& p);
    /**
     * @brief Affectation par égalité
     * @param ilist exemple : {2., 4., 5}
     * @return Point
     */
    Point& operator= (std::initializer_list<double> ilist);

    /**
     * @brief Affectation par égalité
     * @param ilist exemple : {2., 4., 5}
     * @return Point
     */
    Point& operator= (std::initializer_list<int> ilist);

    /**
     * @brief Destructeur
     */
    ~Point ();

    /**
     * @brief Définition de l'indice global du point
     * @param index un entier int
     * @return this Point*
     */
    Point* SetGlobalIndex (int index);

    /**
     * @brief Retourne l'indice global du point dans le maillage auquel il appartient
     * @return index l'indice global
     */
    int GetGlobalIndex () const;

    /**
     * @brief Définition de la localisation du point. Voir l'énumération POINT_LOCATION
     * @param loc tag de type POINT_LOCATION
     * @return rien
     */
    Point* SetLocate (POINT_LOCATION loc);

    /**
     * @brief Retourne le tag du point.
     * @return tag de type POINT_LOCATION
     */
    POINT_LOCATION GetLocate () const;

    /**
     * @brief Nettoie la liste des points voisins de ce point
     * @return this Point*
     */
    void ClearListNeighbours ();

    /**
     * @brief Ajoute le point p à la liste de voisins de ce point
     * @param p point de type Point
     * @return this Point*
     */
    void AddPointNeighbour (Point* p);

    /**
     * @brief Retourne un vecteur de points voisins de ce point
     * @return vecteur de point
     */
    std::vector <Point *> GetListNeighbours ();

    /**
     * @brief Associe un vecteur de points voisins à ce point
     * @param list des Point*
     */
    void SetListNeighbours (std::vector <Point *>& list);

    /**
     * @brief Enlève le point p à la liste de voisins de ce point
     * @param p point de type Point
     */
    void RemoveThisNeighbourPoint (Point* p);

    /**
     * @brief Relie la cellule c à ce point
     * @param c la cellule considérée
     */
    void LinkToCell (Cell* c);

    /**
     * @brief Enlève la connection de ce point à la cellule c
     * @param c la cellule considérée
     */
    void UnlinkToCell (Cell* c);

    /**
     * @brief Retourne la liste des cellules liées à ce point
     * @return Vecteur de Cell*
     */
    std::vector <Cell*> GetLinkedCell () const;

    /**
     * @brief DetachFromAll
     * @return
     */
    Point* DetachFromAll ();

    // Fonction friend
    friend std::ostream & operator<< (std::ostream &out, const Point &p);

    friend double EuclidianDist (const Point& a, const Point& b);
    friend Point operator* (double value, const Point& p);
    friend Point operator* (const Point& p, double value);
    friend Point operator* (const Point& a, const Point& b);

    friend Point operator+ (const Point& a, const Point& b);
    friend Point operator- (const Point& a, const Point& b);

    friend double operator| (const Point& a, const Point& b);

    friend Point operator+ (const Point& p, double value);
    friend Point operator+ (double value, const Point& p);

    friend Point operator- (const Point& p, double value);
    friend Point operator- (double value, const Point& p);

    friend Point operator/ (const Point& p, double value);
    friend Point operator/ (double value, const Point& p);

    friend bool operator< (const Point& a, const Point& b);
    friend bool operator> (const Point& a, const Point& b);

    friend bool operator<= (const Point& a, const Point& b);
    friend bool operator>= (const Point& a, const Point& b);

    friend bool operator== (const Point& a, const Point& b);
    friend bool operator!= (const Point& a, const Point& b);

protected:
    /**
     * @brief Liste des cellules connectées en ce point
     */
    std::vector<Cell *> m_cells;

    /**
    * @brief Liste des points voisins
    */
    std::vector<Point *> m_listNeighbours;

    /**
     * @brief Indice global du point dans le maillage
     */
    int m_globalIndex;

    /**
     * @brief Tag de localisation du point
     */
    POINT_LOCATION m_locate;
};

/**
 * @brief Opérateur d'affichage (surcharge des doubles chevrons)
 * operator << : utiliser pour std::cout << p << std::endl, similaire à printf ("%f, %f, %f", p.x, p.y, p.z).
 * @param p Point
 * @param out flux
 * @return flux
 */
std::ostream& operator<< (std::ostream &out, const Point &p);

/**
 * @brief Calcul de la distance euclidienne entre deux points
 * @param a Point
 * @param b Point
 * @return distance
 */
double EuclidianDist (const Point& a, const Point& b);

/**
 * @brief Multiplication scalaire
 * @param value double
 * @param p le point en question
 * @return Point
 */
Point operator* (double value, const Point& p);

/**
 * @brief Multiplication scalaire
 * @param value double
 * @param p le point en question
 * @return Point
 */
Point operator* (const Point& p, double value);

/**
 * @brief Multiplication intelligente entre deux points
 * @param a le premier point
 * @param b le deuxième point
 * @return Point
 */
Point operator* (const Point& a, const Point& b);

/**
 * @brief Addition de deux points
 * @param a le premier point
 * @param b le deuxième point
 * @return Point
 */
Point operator+ (const Point& a, const Point& b);

/**
 * @brief Soustraction de deux points
 * @param a le premier point
 * @param b le deuxième point
 * @return Point
 */
Point operator- (const Point& a, const Point& b);

/**
 * @brief Produit scalaire entre deux points (vecteurs)
 * @param a le premier point
 * @param b le deuxième point
 * @return double
 */
double operator| (const Point& a, const Point& b);

/**
 * @brief Addition d'un point et d'un double
 * @param p le point en question
 * @param value la valeur à ajouter
 * @return Point
 */
Point operator+ (const Point& p, double value);

/**
 * @brief Addition d'un point et d'un double
 * @param p le point en question
 * @param value la valeur à ajouter
 * @return Point
 */
Point operator+ (double value, const Point& p);

/**
 * @brief Soustraction d'un point et d'un double
 * @param p le point en question
 * @param value la valeur à soustraire
 * @return Point
 */
Point operator- (const Point& p, double value);

/**
 * @brief Soustraction d'un point et d'un double
 * @param p le point en question
 * @param value la valeur à soustraire
 * @return Point
 */
Point operator- (double value, const Point& p);

/**
 * @brief Division à droite d'un point par un double
 * @param p le point en question
 * @param value la valeur en question
 * @return Point
 */
Point operator/ (const Point& p, double value);

/**
 * @brief Division à gauche d'un point par un double
 * @param p le point en question
 * @param value la valeur en question
 * @return Point
 */
Point operator/ (double value, const Point& p);

/**
 * @brief Operateur de comparaison <
 * @param a le premier point
 * @param b le deuxième point
 * @return booléen
 */
bool operator< (const Point& a, const Point& b);

/**
 * @brief Operateur de comparaison >
 * @param a le premier point
 * @param b le deuxième point
 * @return booléen
 */
bool operator> (const Point& a, const Point& b);

/**
 * @brief Operateur de comparaison >=
 * @param a le premier point
 * @param b le deuxième point
 * @return booléen
 */
bool operator<= (const Point& a, const Point& b);

/**
 * @brief Operateur de comparaison >=
 * @param a le premier point
 * @param b le deuxième point
 * @return booléen
 */
bool operator>= (const Point& a, const Point& b);

/**
 * @brief Operateur d'égalité entre deux points (sur les coordonnées)
 * @param a le premier point
 * @param b le deuxième point
 * @return booléen
 */
bool operator== (const Point& a, const Point& b);

/**
 * @brief Operateur de différence entre deux points (sur les coordonnées)
 * @param a le premier point
 * @param b le deuxième point
 * @return booléen
 */
bool operator!= (const Point& a, const Point& b);

#endif // POINT_H
