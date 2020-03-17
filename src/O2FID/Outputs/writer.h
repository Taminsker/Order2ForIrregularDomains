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
     * @brief Constructeur
     */
    Writer (Mesh * mesh);
    /**
     * @brief Destructeur
     */
    ~Writer ();

    /**
     * @brief SetVectorAnalytical
     * @param vector
     */
    void SetVectorAnalytical (Vector * vector);

    /**
     * @brief SetVectorNumerical
     * @param vector
     */
    void SetVectorNumerical (Vector * vector);

    /**
     * @brief SetVectorErrorAbs
     * @param vector
     */
    void SetVectorErrorAbs (Vector * vector);

    /**
     * @brief SetFilename
     * @param filename
     */
    void SetFilename (std::string filename);

    /**
     * @brief WriteBothDomainsOn
     */
    void SetWriteBothDomainsOn ();

    /**
     * @brief WriteBothDomainsOff
     */
    void SetWriteBothDomainsOff ();

    /**
     * @brief SetCurrentIteration
     * @param i
     */
    void SetCurrentIteration (int i);

    /**
     * @brief WriteNow
     */
    void WriteNow ();

protected:

    Mesh * m_mesh;
    /**
     * @brief m_filename
     */
    std::string m_filename;

    /**
     * @brief m_index
     */
    int m_index;

    /**
     * @brief m_domain
     */
    bool m_bothDomain;

    /**
     * @brief m_sol_num
     */
    Vector * m_sol_num;

    /**
     * @brief m_sol_ana
     */
    Vector * m_sol_ana;

    /**
     * @brief m_sol_abs
     */
    Vector * m_error_abs;
};

void WriteInFile (std::ofstream &file, std::string name, Vector * vec);
#endif // WRITER_H
