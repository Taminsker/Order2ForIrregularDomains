/** \file border.h */

#ifndef BORDER_H
#define BORDER_H

#include "../Base/abstractbase.h"
#include "../Data/data.h"
#include "../Mesh/mesh.h"

/**
 * @brief La classe Border sert à lister les points sur la frontière.
 */
class Border : public AbstractBase
{
public:
    /**
     * @brief Construteur.
     */
    Border (Mesh * mesh);

    /**
     * @brief Destructeur.
     */
    ~Border ();

    /**
     * @brief SetFromLevelSet
     */
    void SetLevelSet (double (*phi)(Point));

    /**
     * @brief Création de la liste d'index de points sur le zéro de la fonction
     * @return
     */
    std::vector<int> MakeListOfIndexPoints ();

protected:
    /**
    * @brief Fonction levelSet
    */
    double (*m_phi) (Point);


};

#endif // BORDER_H
