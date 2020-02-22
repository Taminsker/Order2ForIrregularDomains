#ifndef CELL_H
#define CELL_H

#include <vector>

class Point;

/**
 * @brief La classe Cell implémente les informations nécessaires pour définir une cellule.
 * Cette classe met en lien les cellules voisines, les points définissant la cellule et calcule le barycentre.
 */
class Cell
{
public:
    /**
     * @brief Créer une cellule vide
     */
    Cell ();

    /**
     * @brief Détruit la cellule
     */
    ~Cell ();

    /**
     * @brief Configure la cellule au Nord
     * @param cell Cell *
     */
    void SetCellNorth (Cell * cell);

    /**
     * @brief Configure la cellule au Sud
     * @param cell Cell *
     */
    void SetCellSouth (Cell * cell);

    /**
     * @brief Configure la cellule à l'Ouest
     * @param cell Cell *
     */
    void SetCellWest (Cell * cell);

    /**
     * @brief Configure la cellule à l'Est
     * @param cell Cell *
     */
    void SetCellEast (Cell * cell);

    /**
     * @brief Ajout du point p à la cellule
     * @param p Point
     */
    void AddPoint (const Point &p);

    /**
     * @brief Retourne la cellule située au Nord
     * @return Cell *
     */
    Cell * GetCellNorth ();

    /**
     * @brief Retourne la cellule située au Sud
     * @return Cell *
     */
    Cell * GetCellSouth ();

    /**
     * @brief Retourne la cellule située à l'Est
     * @return Cell *
     */
    Cell * GetCellWest ();

    /**
     * @brief Retourne la cellule située à l'Ouest
     * @return Cell *
     */
    Cell * GetCellEast ();

    /**
     * @brief Retourne le barycentre de la cellule
     * @return
     */
    Point GetBarycenter ();

private:
    /**
     * @brief Cellule au Nord
     */
    Cell * m_cell_n;

    /**
     * @brief Cellule au Sud
     */
    Cell * m_cell_s;

    /**
     * @brief Cellule à l'Ouest
     */
    Cell * m_cell_o;

    /**
     * @brief Cellule à l'Est.
     */
    Cell * m_cell_e;

    /**
     * @brief Liste des points de la cellule
     */
    std::vector<Point *> m_points;

    /**
     * @brief Barycentre de la cellule.
     */
    Point * m_bary;
};

#endif // CELL_H
