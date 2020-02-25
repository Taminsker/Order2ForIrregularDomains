/** \file errorsbuilder.h */
#ifndef ERRORSBUILDER_H
#define ERRORSBUILDER_H

#include "../Data/data.h"
#include "../Mesh/mesh.h"

/**
 * @brief La classe ErrorsBuilder sert pour générer les différentes erreurs.
 */
class ErrorsBuilder
{
public:
    /**
     * @brief Constructeur
     */
    ErrorsBuilder (Mesh * mesh);

    /**
     * @brief Destructeur
     */
    ~ErrorsBuilder ();

    /**
     * @brief SetVectorNumerical
     * @param u
     */
    void SetVectorNumerical (Vector * u);

    /**
     * @brief SetVectorAnalytical
     * @param u
     */
    void SetVectorAnalytical (Vector * u);

    /**
     * @brief GetErrorL2
     * @return
     */
    double GetErrorL2 ();

    /**
     * @brief GetErrorLinf
     * @return
     */
    double GetErrorLinf ();

    /**
     * @brief GetErrorRela
     * @return
     */
    double GetErrorRela ();

    /**
     * @brief GetErrorAbs
     * @return
     */
    Vector GetErrorAbs ();

protected:
    /**
     * @brief m_analytical
     */
    Vector * m_analytical;

     /**
      * @brief m_numerical
      */
     Vector * m_numerical;

     Mesh * m_mesh;



};

#endif // ERRORSBUILDER_H
