/** \file impose.h */
#ifndef IMPOSE_H
#define IMPOSE_H

#include "../Base/abstractbase.h"
#include "../Data/data.h"
#include "../Mesh/mesh.h"


/**
 * @brief La classe Impose sert à imposer des conditions de Dirichlet ou de Neumann.
 *
 */
class Impose : public AbstractBase
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
     * @brief MakeNeumann
     * @param sparsematrix
     * @param f
     * @return
     */
    Vector MakeNeumann (Matrix * sparsematrix, double (*f) (Point));

protected:
    /**
     * @brief Vecteur d'index de Points.
     * @return vector
     */
    std::vector <int> m_indexPoints;
};

#endif // IMPOSE_H
