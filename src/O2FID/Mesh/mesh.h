/** @file mesh.h */

#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <iomanip>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../Data/data.h"
#include "../Toolbox/toolbox.h"

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
 * @brief Énumérateur des différentes dimensions possibles
 */
typedef enum {
    /** dimension 1 */
    DIM_1D = 1,
    /** dimension 2 */
    DIM_2D = 2,
    /** dimension 3 */
    DIM_3D = 3
} DIM;

/**
 * @brief Structure embarquée regroupant les indices i, j, k
 */
typedef struct
{
    int i;
    int j;
    int k;
} LocalIndexes;

/** @}*/
/**
 * @brief La classe Mesh : stock un maillage, liste des points et des cellules.
 */
class Mesh
{

public:
    /**
     * @brief Créer un maillage Nx=Ny=Nz=0 et Origin=Extrema=(0, 0, 0).
     */
    Mesh ();

    /**
     * @brief Détruit le maillage.
     */
    ~Mesh ();

    /**
     * @brief Configurer le nombre de points dans la direction x.
     * @param Nx int
     * @return this Mesh*
     */
    Mesh* Set_Nx (int Nx);

    /**
     * @brief Configurer le nombre de points dans la direction y.
     * @param Ny int
     * @return this Mesh*
     */
    Mesh* Set_Ny (int Ny);

    /**
     * @brief Configurer le nombre de points dans la direction z.
     * @param Nz int
     * @return this Mesh*
     */
    Mesh* Set_Nz (int Nz);

    /**
     * @brief Configurer les points définissant le domaine [a, b] x [c, d] x [e, f].
     * @param Origin Point (a, c, e)
     * @param Extrema Point (b, d, f)
     */
    Mesh* SetBounds (Point* Origin, Point* Extrema);

    /**
     * @brief Construit le maillage à partir des données rentrées : Origin, Extrema, Nx, Ny et Nz.
     */
    Mesh* Build ();

    /**
     * @brief Retourne un pointeur vers le point (i, j, k) ie le point situé sur la i-eme ligne, la j-ieme colonne et la k-ieme hauteur
     * @param i indice de la ligne
     * @param j indice de la colonne
     * @param k indice de la hauteur
     * @return Point*
     */
    Point* GetPoint (int i, int j, int k) const;

    /**
     * @brief Retourne un pointeur vers le i-ème point du maillage.
     * @param i int numéro du point
     * @return Point*
     */
    Point* GetPoint (int i) const;

    /**
     * @brief Retourne un pointeur vers le point (i, j, k) ordonné sous forme de parcours en snake
     * @param i indice de la ligne
     * @param j indice de la colonne
     * @param k indice de la hauteur
     * @return Point*
     */
    Point* GetSnakePoint (int i, int j, int k) const;

    /**
     * @brief Retourne un pointeur vers le i-ième point ordonné sous forme de parcours en snake
     * @return Point*
     */
    Point* GetSnakePoint (int idx) const;

    /**
     * @brief GetLocalIndexesOfPoint
     * @param idx
     * @return
     */
    LocalIndexes GetLocalIndexesOfPoint (int idx) const;

    /**
     * @brief GetGlobalIndexOfPoint
     * @param i
     * @param j
     * @param k
     * @return
     */
    int GetGlobalIndexOfPoint (int i, int j = 0, int k = 0) const;


    /**
     * @brief Retourne un pointeur vers la cellule (i, j, k) ie la cellule située sur la i-eme ligne, la j-ieme colonne et la k-ieme hauteur
     * @param i indice de la ligne
     * @param j indice de la colonne
     * @param k indice de la hauteur
     * @return Cell*
     */
    Cell* GetCell (int i, int j, int k) const;

    /**
     * @brief Retourne un pointeur vers la i-ème cellule du maillage.
     * @param i int numéro de la cellule
     * @return Cell*
     */
    Cell* GetCell (int i) const;

    /**
     * @brief GetLocalIndexesOfCell
     * @param idx
     * @return
     */
    LocalIndexes GetLocalIndexesOfCell (int idx) const;

    /**
     * @brief GetGlobalIndexOfCell
     * @param i
     * @param j
     * @param k
     * @return
     */
    int GetGlobalIndexOfCell (int i, int j = 0, int k = 0) const;


    /**
     * @brief Retourne la dimension du maillage (voir l'énumération DIM)
     * @return tag de l'énumération DIM
     */
    DIM GetDimension () const;

    /**
     * @brief Retourne un vecteur contenant des pointeurs Origin et Extrema.
     * @return std::vector<Point*> = {Origin, Extrema}
     */
    std::vector<Point *> GetBounds () const;

    /**
     * @brief Retourne le Nx enregistré.
     * @return Nx int
     */
    int Get_Nx () const;

    /**
     * @brief Retourne le Ny enregistré.
     * @return Ny int
     */
    int Get_Ny () const;

    /**
     * @brief Retourne le Nz enregistré.
     * @return Nz int
     */
    int Get_Nz () const;

    /**
     * @brief Retourne le pas d'espace dans la direction x.
     * @return hx double
     */
    double Get_hx () const;

    /**
     * @brief Retourne le pas d'espace dans la direction y.
     * @return hy double
     */
    double Get_hy () const;

    /**
     * @brief Retourne le pas d'espace dans la direction z.
     * @return hz double
     */
    double Get_hz () const;

    /**
     * @brief Retourne le nombre total de points du maillage.
     * @return N int
     */
    int GetNumberOfTotalPoints () const;

    /**
     * @brief Retourne le nombre de points sur la grille cartésienne.
     * @return N int
     */
    int GetNumberOfCartesianPoints () const;

    /**
     * @brief Retourne le nombre total de cellules du maillage.
     * @return N int
     */
    int GetNumberOfTotalCells () const;

    /**
     * @brief Tague le point (i, j, k)
     * @param i indice de la ligne
     * @param j indice de la colonne
     * @param k indice de la hauteur
     * @param tag un tag de POINT_LOCATION
     * @return this Mesh *
     */
    Mesh* TagPoint (int i, int j, int k, POINT_LOCATION tag);

    /**
     * @brief Tague le index-ieme point
     * @param index indice global du point
     * @param tag un tag de POINT_LOCATION
     * @return this Mesh *
     */
    Mesh* TagPoint (int index, POINT_LOCATION tag);

    /**
     * @brief Ajoute le point au maillage et le tag comme un point de bord (voir POINT_LOCATION)
     * @param a le point considéré
     * @return this Mesh *
     */
    Mesh* AddPointOnBorder (Point a);

    /**
     * @brief Ajoute le ptr point au maillage et le tag comme un point de bord (voir POINT_LOCATION)
     * @param a le point considéré
     * @return this Mesh *
     */
    Mesh* AddPointOnBorder (Point* a);

    /**
     * @brief Ajoute le point au maillage et le tag comme un point de du maillage (voir POINT_LOCATION)
     * @param a le point considéré
     * @param tag voir POINT_LOCATION
     * @return this Mesh *
     */
    Mesh* AddPointOnDomain (Point a, POINT_LOCATION tag = ON_DOMAIN_EXTERN_OMEGA);

    /**
     * @brief Ajoute le ptr point au maillage et le tag comme un point de du maillage (voir POINT_LOCATION)
     * @param a le point considéré
     * @param tag voir POINT_LOCATION
     * @return this Mesh *
     */
    Mesh* AddPointOnDomain (Point* a, POINT_LOCATION tag = ON_DOMAIN_EXTERN_OMEGA);

    /**
     * @brief MakeACellFromListPoints
     * @param list
     * @return
     */
    Mesh* MakeACellFromListPoints (std::vector<Point *> list);
    /**
     * @brief Retourne un vecteur d'indices globaux de points avec le tag 'tag'
     * @param tag voir POINT_LOCATION
     * @return un vecteur d'indices globaux
     */
    std::vector <int> GetListOfIndexPoints (POINT_LOCATION tag = ON_BORDER_OMEGA);

    /**
     * @brief Retourne un vecteur d'indices globaux de cellules avec le tag 'tag'
     * @param tag voir CELL_LOCATION
     * @return un vecteur d'indices globaux
     */
    std::vector <int> GetListOfIndexCells (CELL_LOCATION tag = IN_DOMAIN_INTERN_OMEGA);

    /**
     * @brief Retourne le nombre total d'information sur les cellules disponibles (fonction utilisée pour écrire un fichier VTK, partie CELLS)
     * @return N un entier
     */
    int GetNumberOfInfosCells () const;

    /**
     * @brief Élimine du maillage le point index
     * @param index indice global du point
     * @return this Mesh*
     */
    Mesh* RemoveThisPoint (int index);

    /**
     * @brief Élimine du maillage le point (i, j, k)
     * @param i indice de la ligne
     * @param j indice de la colonne
     * @param k indice de la hauteur
     * @return this Mesh*
     */
    Mesh* RemoveThisPoint (int i, int j, int k);

    /**
     * @brief Élimine du maillage la liste de points donné par indices
     * @param listindex list de indices
     * @return this Mesh*
     */
    Mesh* RemoveThesePoints (std::vector <int> listindex);

    /**
     * @brief Élimine tous les points hors grille cartesienne
     * @return this Mesh*
     */
    Mesh* RemoveAllNotCartesianPoints ();

    /**
     * @brief Affiche à l'écran les infos et la liste de points.
     * @return this
     */
    Mesh* Print ();

    /**
     * @brief Met le vecteur à zéro sur les indices de points à l'extèrieur du domaine
     * @param vec pointeur vers un Vector
     */
    void MakeZeroOnExternOmegaInVector (Vector * vec);


protected:

    /**
     * @brief Vecteur de Point.
     */
    std::vector<Point*> m_points;

    /**
     * @brief Vecteur de Cell.
     */
    std::vector<Cell*> m_cells;

    /**
     * @brief Point d'origine du maillage.
     */
    Point* m_origin;

    /**
     * @brief Point d'extrema du maillage.
     */
    Point* m_extrema;

    /**
     * @brief Pas d'espace dans la direction x.
     */
    double m_hx;

    /**
     * @brief Pas d'espace dans la direction y.
     */
    double m_hy;

    /**
     * @brief Pas d'espace dans la direction z.
     */
    double m_hz;

    /**
     * @brief Nombre de points dans la direction x.
     */
    int m_Nx;

    /**
     * @brief Nombre de points dans la direction y.
     */
    int m_Ny;

    /**
     * @brief Nombre de points dans la direction z.
     */
    int m_Nz;

    /**
     * @brief Retourne la position du point dans le vecteur de Point de numéro (i, j, k) ie son indice global.
     * @param i int ordre dans la direction x
     * @param j int ordre dans la direction y
     * @param k int ordre dans la direction z
     * @return index int indice global
     */
    int IndexPoints (int i, int j, int k) const;

    /**
     * @brief Retourne la position de la cellule dans le vecteur de Cell de numéro (i, j, k) ie son indice global.
     * @param i int ordre dans la direction x
     * @param j int ordre dans la direction y
     * @param k int ordre dans la direction z
     * @return index int indice global
     */
    int IndexCells (int i, int j, int k) const;
};

#endif // MESH_H
