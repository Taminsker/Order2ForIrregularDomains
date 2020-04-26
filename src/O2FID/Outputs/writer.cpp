#include "writer.h"

Writer::Writer (Mesh * mesh) :
    m_mesh (mesh),
    m_filename ("no_filename_selected.vtk"),
    m_index (0),
    m_bothDomain (true),
    m_writeDat (false),
    m_sol_num (nullptr),
    m_sol_ana (nullptr),
    m_error_abs (nullptr),
    m_phi_value (nullptr),
    m_normals (nullptr),
    m_Wnew (nullptr),
    m_Wold (nullptr)
{}

Writer::~Writer ()
{}
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
    m_base_filename = filename;
    return;
}

void Writer::SetVectorPhi (Vector *vector)
{
    m_phi_value = vector;
    return;
}
void Writer::SetVectorNormals (std::vector<Point *> *vector)
{
    m_normals = vector;
    return;
}

void Writer::SetVectorW_new (std::vector<Point *>* vector)
{
    delete m_Wnew;
    m_Wnew = new std::vector<Point *>(*vector);
    return;
}

void Writer::SetVectorW_old (std::vector<Point *>* vector)
{
    delete m_Wold;
    m_Wold = new std::vector<Point *>(*vector);
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
    std::cout << "# The writer is launched." << std::endl;

    if (m_filename == "")
        m_filename = "no_filename_selected";

    std::cout << INDENT << "Base filename is " << m_base_filename << std::endl;
    std::cout << INDENT << "Option : bothDomain is " << (m_bothDomain ? "on": "off") << std::endl;

    // Ajout du numéro d'itération en temps et de l'extension du fichier ".vtk"
    m_filename = m_base_filename + std::string ("_") + std::to_string (m_index);

    if (m_writeDat)
        WriteDAT ();
    WriteVTK ();

    return;
}

void Writer::WriteDAT ()
{    
    std::ofstream file;
    std::string filename_copy;

    std::vector<int> indexes (size_t (m_mesh->GetNumberOfTotalPoints ()));

    if (m_bothDomain)
    {
        for (int i = 0; i < m_mesh->GetNumberOfTotalPoints (); ++i)
            indexes.at (size_t (i)) = i;
    }
    else
    {
        indexes = m_mesh->GetListOfIndexPoints ();
    }

    int numPoints = int (indexes.size ());

    // POINTS
    filename_copy = m_filename + std::string ("_points.dat");

    std::cout << INDENT << "(DAT) Write : " << filename_copy << std::endl;

    file.open (filename_copy);

    if (file.is_open ())
    {
        file << "# X " << SPACE << " Y " << SPACE << " Z " << std::endl;

        // Impression de la liste de points du maillage
        for (int i : indexes)
            file << *m_mesh->GetPoint (i) << std::endl;

        file.close ();
    }

    filename_copy = m_filename + std::string ("_location.dat");

    std::cout << INDENT << "(DAT) Write : " << filename_copy << std::endl;

    file.open (filename_copy);


    // Vecteur de localisation (Domaine externe, frontière et domaine interne).
    file << "#X" << SPACE << "Y" << SPACE << "Z" << SPACE << "LOC" << std::endl;
    for (int i : indexes)
        file << *m_mesh->GetPoint (i) << SPACE << m_mesh->GetPoint (i)->GetLocate () << std::endl;

    file.close ();


    // Impression du vecteur de solution numérique si disponible
    if (m_sol_num != nullptr && m_sol_num->rows () >= numPoints)
    {
        filename_copy = m_filename + std::string ("_sol_num.dat");

        std::cout << INDENT << "(DAT) Write : " << filename_copy << std::endl;

        file.open (filename_copy);

        file << "#X" << SPACE << "Y" << SPACE << "Z" << SPACE << "SOL_NUM" << std::endl;
        for (int i : indexes)
            file << *m_mesh->GetPoint (i) << SPACE << m_sol_num->operator() (i) << std::endl;

        file.close ();
    }

    // Impression du vecteur de solution analytique si disponible
    if (m_sol_ana != nullptr && m_sol_ana->rows () >= numPoints)
    {
        filename_copy = m_filename + std::string ("_sol_ana.dat");

        std::cout << INDENT << "(DAT) Write : " << filename_copy << std::endl;

        file.open (filename_copy);

        file << "#X" << SPACE << "Y" << SPACE << "Z" << SPACE << "SOL_ANA" << std::endl;
        for (int i : indexes)
            file << *m_mesh->GetPoint (i) << SPACE << m_sol_ana->operator() (i) << std::endl;

        file.close ();
    }


    // Impression du vecteur d'erreurs en valeurs absolues si disponibles
    if (m_error_abs != nullptr && m_error_abs->rows () >= numPoints)
    {
        filename_copy = m_filename + std::string ("_error_abs.dat");

        std::cout << INDENT << "(DAT) Write : " << filename_copy << std::endl;

        file.open (filename_copy);

        file << "#X" << SPACE << "Y" << SPACE << "Z" << SPACE << "ERR_ABS" << std::endl;
        for (int i : indexes)
            file << *m_mesh->GetPoint (i) << SPACE << m_error_abs->operator() (i) << std::endl;

        file.close ();
    }

    std::cout << std::endl;

    return;
}

void Writer::WriteVTK ()
{
    std::ofstream file;

    std::string filename_copy = m_filename;
    filename_copy += std::string (".vtk");

    std::cout << INDENT << "(VTK) Write : " << filename_copy << std::endl;

    file.open (filename_copy);

    int numPoints = m_mesh->GetNumberOfTotalPoints ();

    file << "# vtk DataFile Version 2.0" << std::endl;
    file << filename_copy << ", O2FID output." << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << numPoints << " double" << std::endl;

    // Impression de la liste de points du maillage
    for (int i = 0; i < numPoints; ++i)
        file << *m_mesh->GetPoint (i) << std::endl;

    file << std::endl;

    std::cout << INDENT << "(VTK) Points have been written." << std::endl;

    // Partie concernant l'impression des cellules

    int numCells = m_mesh->GetNumberOfTotalCells ();
    int numInfosCells = m_mesh->GetNumberOfInfosCells ();

    if (m_bothDomain) // Si on a précisé qu'on voulait les deux domaines
    {
        file << "CELLS " << numCells << " " << numInfosCells << std::endl;

        // Impression de la liste des cellules
        for (int i = 0; i < numCells; ++i)
            file << *m_mesh->GetCell (i);
        file << std::endl;

        std::cout << INDENT << "(VTK) Cells have been written." << std::endl;

        file << "CELL_TYPES " << numCells << std::endl;

        // Impression de la liste des types de cellules (voir documentation VTK pour plus d'information)
        for (int i = 0; i < numCells; ++i)
            file << m_mesh->GetCell (i)->GetType () << std::endl;
        file << std::endl;

        std::cout << INDENT << "(VTK) Cell type was written," << std::endl;

    }
    else // On a précisé qu'on ne voulait que le domaine interne
    {
        std::vector <int> indexCells = m_mesh->GetListOfIndexCells ();
        numCells = int (indexCells.size ());

        // Récupère le nombre d'infos des cellules

        numInfosCells = 0;
        for (int i = 0; i < numCells; ++i)
        {
            int index = indexCells.at (size_t (i));
            numInfosCells += m_mesh->GetCell (index)->GetNumberOfInfos ();
        }

        // Écrire les cellules taguées
        file << "CELLS " << numCells << " " << numInfosCells << std::endl;
        for (int i = 0; i < numCells; ++i)
        {
            int index = indexCells.at (size_t (i));
            file << *m_mesh->GetCell (index);
        }
        file << std::endl;

        std::cout << INDENT << "(VTK) Cells have been written." << std::endl;

        // Écrire le type des cellules taguées
        file << "CELL_TYPES " << numCells << std::endl;

        for (int i = 0; i < numCells; ++i)
        {
            int index = indexCells.at (size_t (i));
            file << m_mesh->GetCell (index)->GetType () << std::endl;
        }
        file << std::endl;

        std::cout << INDENT << "(VTK) Cell type was written." << std::endl;
    }

    // Impression des vecteurs de données sur les points du maillage dans des vecteurs scalaires.
    file << "POINT_DATA " << numPoints << std::endl;

    // Vecteur de localisation (Domaine externe, frontière et domaine interne).
    Vector loc (numPoints);
    loc.setOnes ();

    for (int i = 0; i < numPoints; ++i)
        loc (i) = m_mesh->GetPoint (i)->GetLocate ();

    // Impression du vecteur de localisation
    WriteInFile (file, "Location", &loc);

    // Impression du vecteur de solution numérique si disponible
    if (m_sol_num != nullptr && m_sol_num->rows () == numPoints)
        WriteInFile (file, "Sol_num", m_sol_num);

    // Impression du vecteur de solution analytique si disponible
    if (m_sol_ana != nullptr && m_sol_ana->rows () == numPoints)
        WriteInFile (file, "Sol_ana", m_sol_ana);

    // Impression du vecteur d'erreurs en valeurs absolues si disponibles
    if (m_error_abs != nullptr && m_error_abs->rows () == numPoints)
        WriteInFile (file, "Error_abs", m_error_abs);

    if (m_phi_value != nullptr && m_phi_value->rows () == numPoints)
        WriteInFile (file, "Phi_Value", m_phi_value);

    if (m_normals != nullptr && int(m_normals->size ()) == numPoints)
    {
        file << "NORMALS Normal double" << std::endl;
        for (size_t i = 0; i < m_normals->size (); ++i)
            file << *m_normals->at (i) << std::endl;

        file << std::endl;

        std::cout << INDENT << "(VTK) POINT_DATA Normal was written." << std::endl;

    }

    if (m_Wnew != nullptr && int(m_Wnew->size ()) == numPoints)
    {
        file << "VECTORS W_Current double" << std::endl;
        for (size_t i = 0; i < m_Wnew->size (); ++i)
            file << *m_Wnew->at (i) << std::endl;

        file << std::endl;

        std::cout << INDENT << "(VTK) POINT_DATA W_Current was written." << std::endl;
    }

    if (m_Wold != nullptr && int(m_Wold->size ()) == numPoints)
    {
        file << "VECTORS W_Old double" << std::endl;
        for (size_t i = 0; i < m_Wold->size (); ++i)
            file << *m_Wold->at (i) << std::endl;

        file << std::endl;

        std::cout << INDENT << "(VTK) POINT_DATA W_Old was written." << std::endl;
    }

    file.close ();

    std::cout << std::endl;

    return;
}

void WriteInFile (std::ofstream &file, std::string name, Vector * vec)
{
    file << "SCALARS " << name << " double" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;

    // Impression du vecteur vec
    for (int i = 0; i < vec->rows (); ++i)
        file << vec->operator() (i) << std::endl;

    file << std::endl;

    std::cout << INDENT << "(VTK) POINT_DATA " << name << " was written." << std::endl;

    return;
}
