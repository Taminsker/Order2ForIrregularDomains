#include "writer.h"

Writer::Writer (Mesh * mesh) :
    m_mesh (mesh),
    m_filename ("no_filename_selected.vtk"),
    m_index (0),
    m_bothDomain (true),
    m_sol_num (nullptr),
    m_sol_ana (nullptr),
    m_error_abs (nullptr)
{

}

Writer::~Writer ()
{

}
void Writer::SetVectorAnalytical (Vector *vector)
{
    m_sol_ana = vector;
    return;
}

void Writer::SetVectorNumerical (Vector *vector)
{
    m_sol_num = vector;
    return;
}

void Writer::SetVectorErrorAbs (Vector *vector)
{
    m_error_abs = vector;
    return;
}

void Writer::SetFilename (std::string filename)
{
    m_filename = filename;
    return;
}

void Writer::SetWriteBothDomainsOn ()
{
    m_bothDomain = true;
    return;
}

void Writer::SetWriteBothDomainsOff ()
{
    m_bothDomain = false;
    return;
}

void Writer::SetCurrentIteration (int i)
{
    m_index = i;
    return;
}

void Writer::WriteNow ()
{
    std::ofstream file;
    if (m_filename == "")
        m_filename = "no_filename_selected";

    m_filename += std::string ("_") + std::to_string (m_index) + std::string (".vtk");
    file.open (m_filename);

    int numPoints = m_mesh->GetNumberOfTotalPoints ();

    file << "# vtk DataFile Version 2.0" << std::endl;
    file << m_filename << ", O2FID output." << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << numPoints << " double" << std::endl;

    for (int i = 0; i < numPoints; ++i)
        file << m_mesh->GetPoint (i) << std::endl;

    file << std::endl;

    int numCells = m_mesh->GetNumberOfTotalCells ();
    int numInfosCells = m_mesh->GetNumberOfInfosCells ();

    if (m_bothDomain)
    {
        file << "CELLS " << numCells << " " << numInfosCells << std::endl;

        for (int i = 0; i < numCells; ++i)
            file << m_mesh->GetCell (i);
        file << std::endl;

        file << "CELL_TYPES " << numCells << std::endl;

        for (int i = 0; i < numCells; ++i)
            file << m_mesh->GetCell (i).GetType () << std::endl;
        file << std::endl;
    }
    else
    {
        std::vector <int> indexCells = m_mesh->GetListOfIndexCells ();
        numCells = int (indexCells.size ());

        // Récupère le nombre d'infos des cellules

        numInfosCells = 0;
        for (int i = 0; i < numCells; ++i)
        {
            int index = indexCells.at (size_t (i));
            numInfosCells += m_mesh->GetCell (index).GetNumberOfInfos ();
        }

        // Écrire les cellules taguées
        file << "CELLS " << numCells << " " << numInfosCells << std::endl;
        for (int i = 0; i < numCells; ++i)
        {
            int index = indexCells.at (size_t (i));
            file << m_mesh->GetCell (index);
        }
        file << std::endl;

        // Écrire le type des cellules taguées
        file << "CELL_TYPES " << numCells << std::endl;

        for (int i = 0; i < numCells; ++i)
        {
            int index = indexCells.at (size_t (i));
            file << m_mesh->GetCell (index).GetType () << std::endl;
        }
        file << std::endl;
    }

    file << "POINT_DATA " << numPoints << std::endl;

    Vector loc (numPoints);
    loc.setOnes ();

    for (int i = 0; i < numPoints; ++i)
        loc (i) = m_mesh->GetPoint (i).GetLocate ();

    WriteInFile (file, "Location", &loc);

    if (m_sol_num != nullptr && m_sol_num->rows () == numPoints)
        WriteInFile (file, "Sol_num", m_sol_num);

    if (m_sol_ana != nullptr && m_sol_ana->rows () == numPoints)
        WriteInFile (file, "Sol_ana", m_sol_ana);

    if (m_error_abs != nullptr && m_error_abs->rows () == numPoints)
        WriteInFile (file, "Error_abs", m_error_abs);

    return;
}

void WriteInFile (std::ofstream &file, std::string name, Vector * vec)
{
    file << "SCALARS " << name << " double" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;

    for (int i = 0; i < vec->rows (); ++i)
        file << vec->operator() (i) << std::endl;

    file << std::endl;
    return;
}
