/** \file writer.h */
#ifndef WRITER_H
#define WRITER_H

#include <string>
#include <fstream>

#include "../Data/data.h"
#include "../Mesh/mesh.h"

/**
 * @brief La classe Writer sert pour créer les fichiers vtk et data de sortie.
 */
class Writer
{
public:
    /**
     * @brief Constructeur par initilisation, contruction autour du maillage mesh.
     * @param mesh pointeur vers un maillage de type Mesh.
     */
    Writer (Mesh * mesh);

    /**
     * @brief Destructeur simple.
     */
    ~Writer ();

    /**
     * @brief Définition du vecteur qui contient la solution analytique.
     * @param vector qui est un Vector (voir le fichier datatypedefintions.h).
     * @return rien.
     */
    void SetVectorAnalytical (Vector * vector);

    /**
     * @brief Définition du vecteur qui contient la solution numérique
     * @param vector qui est un Vector (voir le ficher datatypedefintions.h)
     * @return rien
     */
    void SetVectorNumerical (Vector * vector);

    /**
     * @brief Définition du vecteur qui contient l'erreur en valeur absolue sur le maillage entre la solution numérique et la solution analytique
     * @param vector qui est un Vector (voir le ficher datatypedefintions.h)
     * @return rien
     */
    void SetVectorErrorAbs (Vector * vector);

    /**
     * @brief Définition du nom de fichier de sortie utilisé pour écrire le fichier VTK
     * @param filename qui est un std::string
     * @return rien
     */
    void SetFilename (std::string filename);

    void SetVectorPhi (Vector* vector);

    void SetVectorNormals(std::vector<Point*>* vector);

    void SetVectorW_new(std::vector<Point*>* vector);

    void SetVectorW_old(std::vector<Point*>* vector);
    /**
     * @brief On dit que le fichier écrit en sortie comportera le domaine entier sur lequel on travaille
     * @return rien
     */
    void SetWriteBothDomainsOn ();

    /**
     * @brief On dit que le fichier écrit en sortie comportera uniquement le domaine interne lequel on travaille, ie le domaine défini par la fonction levelset (ou phi).
     * @return rien
     */
    void SetWriteBothDomainsOff ();

    /**
     * @brief Les fichiers ".dat" seront écrit
     */
    void SetWriterDatFileOn();

    /**
     * @brief Les fichiers ".dat" ne seront pas écrit
     */
    void SetWriterDatFileOff();


    /**
     * @brief Définition du numéro de l'itération courante (apparaitra dans le nom de fichier en sortie).
     * @param i un entier int
     * @return rien
     */
    void SetCurrentIteration (int i);

    /**
     * @brief Fonction centrale de cette classe. L'écriture du ou des fichiers c'est maintenant.
     * @return rien
     */
    void WriteNow ();

protected:

    /**
     * @brief Pointeur vers un objet de type Mesh.
     */
    Mesh * m_mesh;

    /**
     * @brief Nom du fichier de sortie VTK (sans le ".vtk").
     */
    std::string m_filename;

    /**
     * @brief Nom du fichier de sortie VTK (sans le ".vtk").
     */
    std::string m_base_filename;

    /**
     * @brief L'indice de l'tération de temps courante
     */
    int m_index;

    /**
     * @brief Booléen qui dit si ou non on a choisi d'écrire le domaine en entier ou pas.
     */
    bool m_bothDomain;

    /**
     * @brief Booléen qui dit si ou non on a choisi d'écrire les fichiers ".dat" ou pas
     */
    bool m_writeDat;

    /**
     * @brief Pointeur vers un vecteur de solution numérique.
     */
    Vector * m_sol_num;

    /**
     * @brief Pointeur vers un vecteur de solution analytique.
     */
    Vector * m_sol_ana;

    /**
     * @brief Pointeur vers un vecteur d'erreurs absolues.
     */
    Vector * m_error_abs;

    /**
     * @brief Pointeur vers un vecteur de valeurs de phi
     */
    Vector* m_phi_value;

    /**
     * @brief Pointeur vers un vecteur de normals
     */
    std::vector<Point*>* m_normals;

    /**
     * @brief Pointeur vers un vecteur de champs W au temps courant
     */
    std::vector<Point*>* m_Wnew;

    /**
     * @brief Pointeur vers un vecteur de champs W au temps precédent;
     */
    std::vector<Point*>* m_Wold;


    /**
     * @brief WriteVTK
     */
    void WriteVTK ();

    /**
     * @brief WriteDAT
     */
    void WriteDAT ();
};

/**
 * @brief Fonction qui sert pour écrire les POINT_DATA dans un fichier VTK, voir la doc vtk pour plus d'information sur la structure d'un fichier vtk
 * @param file est un flux ofstream, un flux de sortie vers un fichier
 * @param name est le nom de SCALAR à ajouter (et écrire) dans le fichier VTK
 * @param vec pointeur vers le vecteur à écrire dans le fichier
 * @return rien
 */
void WriteInFile (std::ofstream &file, std::string name, Vector * vec);

#endif // WRITER_H
