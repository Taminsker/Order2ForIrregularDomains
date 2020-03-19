#ifndef CELL_H
#define CELL_H

#include <fstream>
#include <vector>

/**
  * @brief Énumerateur pour la localisation des cellules.
  */
typedef enum {
    IN_DOMAIN_INTERN_OMEGA,
    IN_DOMAIN_EXTERN_OMEGA,
    IN_MIX_DOMAIN
} CELL_LOCATION;

/**
 * @brief Rappel de l'existence de la classe Point (le header n'est pas inclu directement car il y avait un risque d'inclusion réciproque au moment où je codais cette classe).
 */
class Point;

/**
 * @brief La classe Cell implémente les informations nécessaires pour définir une cellule.
 * Cette classe permet nottemment de lier des points au sein d'un même conteneur. Permettant en particulier d'avoir les connectivités (partielles du moins)
 */
class Cell
{
public:
    /**
     * @brief Constructeur par défaut
     */
    Cell ();

    /**
      * @brief Destructeur
      */
    ~Cell ();

    /**
     * @brief Lier un nouveau point avec cette cellule (un pointeur vers ce point est ajouter au sein de la cellule)
     * @param p est un pointeur vers un point de type Point
     * @return rien
     */
    void AddPoint (Point * p);

    /**
     * @brief Désolidariser un point de cette cellule. En pratique supprime juste le pointeur vers le dit point.
     * @param p est un pointeur vers un point de type Point
     * @return rien
     */
    void RemovePoint (Point * p);

    /**
     * @brief Retourne le tag de localisation de la cellule en fonction du tag des points la définissant
     * @return tag de type CELL_LOCATION
     */
    CELL_LOCATION GetLocate () const;

    /**
     * @brief Retourne le type VTK de la cellule (voir documentation VTK_CELL_TYPE pour plus d'informations)
     * @return Type de cellule VTK
     */
    int GetType () const;

    /**
     * @brief Fonction interne pour l'écriture d'un fichier VTK (permet de sommer plus rapidement le nombre de labels qui vont être écris dans le chamsp CELLS du fichier de sortie)
     * @return Nombre de labels de sortie : nbPoints + 1 (label de comptage)
     */
    int GetNumberOfInfos () const;

    /**
     * @brief operator << utilisation pour afficher les numéros globaux des points définissant la cellule
     * @param out flux std::ostream
     * @param c une cellule
     * @return out le flux std::ostream passée en parmamètre
     */
    friend std::ostream & operator<< (std::ostream &out, const Cell &c);

protected:
    /**
    * @brief Vecteur de pointeurs vers les points définissant la cellule
    */
    std::vector <Point *> m_points;

};

/**
 * @brief operator << utilisation pour afficher les numéros globaux des points définissant la cellule
 * @param out flux std::ostream
 * @param c une cellule
 * @return out le flux std::ostream passée en parmamètre
 */
std::ostream & operator<< (std::ostream &out, const Cell &c);


#endif // CELL_H
