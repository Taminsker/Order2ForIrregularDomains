/** @brief matrixbuilder.h */
#ifndef MATRIXBUILDER_H
#define MATRIXBUILDER_H

#include "../Data/data.h"
#include "../Mesh/mesh.h"

/**
 * @brief La classe MatrixBuilder sert à constuire les matrices 1D, 2D, 3D des problèmes de l'équation de la chaleur et de l'équation des Poissons
 */
class MatrixBuilder
{
public:
    /**
     * @brief Constructeur
     */
    MatrixBuilder (Mesh * mesh);

    /**
     * @brief Destructeur
     */
    ~MatrixBuilder ();

    /**
     * @brief SetDelta
     * @param dt
     */
    void SetDeltaT (double dt);

    /**
     * @brief Set1D
     */
    void Set1D ();

    /**
     * @brief Set2D
     */
    void Set2D ();

    /**
     * @brief Set3D
     */
    void Set3D ();

    /**
     * @brief MakeLaplaceEquation
     * @return
     */
    Matrix MakeLaplaceEquation ();

    /**
     * @brief MakeHeatEquation
     * @param coeff
     * @return
     */
    Matrix MakeHeatEquation (double coeff);


protected:

    Mesh * m_mesh;
    /**
     * @brief dim
     */
    int dim = 1;

    /**
     * @brief dt
     */
    double dt;
};

#endif // MATRIXBUILDER_H
