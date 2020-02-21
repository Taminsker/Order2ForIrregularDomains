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

/**
 * @brief Abréviation du type Matrix
 */
typedef Eigen::SparseMatrix<double> Matrix;

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
     */
    void Set_Nx (int Nx);

    /**
     * @brief Configurer le nombre de points dans la direction y.
     * @param Ny int
     */
    void Set_Ny (int Ny);

    /**
     * @brief Configurer le nombre de points dans la direction z.
     * @param Nz int
     */
    void Set_Nz (int Nz);

    /**
     * @brief Configurer les points définissant le domaine [a, b] x [c, d] x [e, f].
     * @param Origin Point : (a, c, e)
     * @param Extrema Point : (b, d, f)
     */
    void SetBounds (Point Origin, Point Extrema);

    /**
     * @brief Construit le maillage à partir des données rentrées : Origin, Extrema, Nx, Ny et Nz.
     */
    void Build ();

    /**
     * @brief Retourne le i-ème point du maillage.
     * @param i int : numéro du point
     * @return Point
     */
    Point operator() (int i) const;

    /**
     * @brief Retourne un vecteur contenant Origin et Extrema.
     * @return std::vector<Point> = {Origin, Extrema}
     */
    std::vector<Point> GetBounds () const;

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
     * @brief Affiche à l'écran les infos et la liste de points.
     */
    void Print ();

protected:

    std::vector<Point> m_points; /**< Vecteur de Point. **/
    Point m_origin; /**< Point d'origine du maillage. **/
    Point m_extrema; /**< Point d'extrema du maillage. **/

    double m_hx; /**< Pas d'espace dans la direction x. **/
    double m_hy; /**< Pas d'espace dans la direction y. **/
    double m_hz; /**< Pas d'espace dans la direction z. **/

    int m_Nx; /**< Nombre de points dans la direction x. **/
    int m_Ny; /**< Nombre de points dans la direction y. **/
    int m_Nz; /**< Nombre de points dans la direction z. **/

    /**
     * @brief Retourne la position du point dans le vecteur de Point de numéro (i, j, k).
     * @param i int ordre dans la direction x
     * @param j int ordre dans la direction y
     * @param k int ordre dans la direction z
     * @return index int index global
     */
    int Index (int i, int j, int k);
};

#endif // MESH_H
