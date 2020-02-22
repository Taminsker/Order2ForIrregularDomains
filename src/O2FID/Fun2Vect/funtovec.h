/** \file funtovec.h */
#ifndef FUNTOVEC_H
#define FUNTOVEC_H

#include "../Base/abstractbase.h"
#include "../Data/data.h"
#include "../Mesh/mesh.h"

/**
 * @brief La classe SourceTermBuilder sert pour créer le terme source à partir d'un maillage et d'une fonction
 */
class FunToVec : public AbstractBase
{
public:
    /**
     * @brief Constructeur
     */
    FunToVec (Mesh * mesh);

    /**
     * @brief Destructeur
     */
    ~FunToVec ();

    /**
     * @brief Make
     * @return
     */
    Vector Make (double (*f) (Point));

    /**
     * @brief Make
     * @return
     */
    Vector Make (double (*f) (Point, double));
    /**
     * @brief Make
     * @param value
     * @return
     */
    Vector Make (double value);
};

#endif // FUNTOVEC_H
