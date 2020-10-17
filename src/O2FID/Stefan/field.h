/** @file field.h */

#ifndef FIELD_H
#define FIELD_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"
#include "../Data/datatypedefinitions.h"
#include <vector>

/*!
 *  \addtogroup Outils
 *  @{
 */

/**
 * @brief Classe utilisée lors du calcul de W sur le domaine, elle regroupe le vecteur des W, le vecteur des vecteurs normaux, le vecteur des gradients de la solution numérique (température) et le vecteur des gradient de phi.
 */
class Field
{
public:
    /**
     * @brief Contructeur
     */
    Field();

    /**
     * @brief Constructeur par copie.
     * @param f l'objet Field à copier.
     */
    Field(const Field& f);

    /**
      * @brief Destructeur
      */
    ~Field();

    /**
     * @brief Surcharge de l'opérateur d'affectation (transfert des données privées).
     * @param f l'objet Field à copier.
     * @return Une référence vers un nouvelle objet Field.
     */
    Field& operator=(const Field& f);


    /**
     * @brief Le champs vectoriel W sur le domaine.
     */
    std::vector<Point*> W;

    /**
     * @brief Les vecteurs normaux à l'interface sur le domaine.
     */
    std::vector<Point*> Normals;

    /**
     * @brief Les gradients de la solution numérique sur le domaine.
     */
    std::vector<Point*> GradTemperature;

    /**
     * @brief Les gradients de la fonction phi (distance) sur le domaine.
     */
    std::vector<Point*> GradPhi;
};
/** @} */

#endif // FIELD_H
