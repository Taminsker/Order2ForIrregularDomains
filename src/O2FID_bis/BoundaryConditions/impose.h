/** \file impose.h */
#ifndef IMPOSE_H
#define IMPOSE_H

#include "../Data/data.h"
#include "../Mesh/mesh.h"


/**
 * @brief La classe Impose sert à imposer des conditions de Dirichlet ou de Neumann.
 *
 */
class Impose
{
public:
    /**
     * @brief Constructeur
     */
    Impose (Mesh * mesh);

    /**
     * @brief Destructeur
     */
    ~Impose ();

    /**
     * @brief Configure la liste des points situés sur le bord Gamma.
     * @param listIndex std::vector <int>
     */
    void SetIndexPoints (std::vector <int> listIndex);

    /**
     * @brief MakeDirichlet
     * @param sparsematrix
     * @param f
     * @return
     */
    Vector MakeDirichlet (Matrix * sparsematrix, double (*f) (Point));

    /**
     * @brief MakeDirichlet
     * @param sparsematrix
     * @param t
     * @return
     */
    Vector MakeDirichlet (Matrix * sparsematrix, double (*f) (Point, double), double t);

    /**
     * @brief MakeNeumann
     * @param sparsematrix
     * @param f
     * @return
     */
    Vector MakeNeumann (Matrix * sparsematrix, double (*f) (Point));

    /**
     * @brief MakeNeumann
     * @param sparsematrix
     * @param t
     * @return
     */
    Vector MakeNeumann (Matrix * sparsematrix, double (*f) (Point, double), double t);


protected:
    /**
     * @brief Vecteur d'index de Points.
     * @return vector
     */
    std::vector <int> m_indexPoints;

    Mesh * m_mesh;
};

#endif // IMPOSE_H
