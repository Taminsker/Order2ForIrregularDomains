/** \file abstractbase.h */
#ifndef ABSTRACTBASE_H
#define ABSTRACTBASE_H

#include "../Data/data.h"
#include "../Mesh/mesh.h"

/**
 * @brief La classe AbstractBase sert de brique de base pour la constrution des différents outils.
 */
class AbstractBase
{
public:
    /**
     * @brief Constructeur
     */
    AbstractBase (Mesh * mesh);

    /**
     * @brief Destructeur
     */
    ~AbstractBase ();

protected:
    /**
     * @brief mesh
     */
    Mesh * m_mesh;
};

#endif // ABSTRACTBASE_H
