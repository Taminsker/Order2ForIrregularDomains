/** @file mesh.cpp */

#include "mesh.h"

Mesh::Mesh () :
    m_origin (new Point ()), m_extrema (new Point ()),
    m_hx (0.), m_hy (0.), m_hz (0.),
    m_Nx (1), m_Ny (1), m_Nz (1)
{}

Mesh::~Mesh ()
{
    for (auto c : m_cells)
        delete c;
    for (auto p : m_points)
        delete p;

    m_cells.clear();
    m_points.clear ();

    delete m_origin;
    delete m_extrema;
}

Mesh* Mesh::Set_Nx (int Nx)
{
    m_Nx = std::abs (Nx);
    return this;
}

Mesh* Mesh::Set_Ny (int Ny)
{
    m_Ny = std::abs (Ny);
    return this;
}

Mesh* Mesh::Set_Nz (int Nz)
{
    m_Nz = std::abs (Nz);
    return this;
}

Mesh* Mesh::SetBounds (Point* Origin, Point* Extrema)
{
    delete m_origin;
    delete m_extrema;

    m_origin = Origin;
    m_extrema = Extrema;
    return this;
}

Mesh* Mesh::Build ()
{
    if (m_Nx == 0 || m_Ny == 0 || m_Nz == 0) // Pas de points
    {
        std::cout << "ERROR Mesh:: triplet for {Nx, Ny, Nz} "
                  << m_Nx << " " << m_Ny << " " << m_Nz << std::endl;
        exit(0);
    }

    if (*m_origin == *m_extrema) // /!\ domaine réduit à un point
    {
        std::cout << "ERROR Mesh:: vector bounds for "
                  << "[" << m_origin->x << ", " << m_extrema->x << "] x "
                  << "[" << m_origin->y << ", " << m_extrema->y << "] x "
                  << "[" << m_origin->z << ", " << m_extrema->z << "]"
                  << std::endl;
        m_Nx = m_Ny = m_Nz = 0;
        exit(0);
    }

    m_hx = (m_extrema->x - m_origin->x) / double (std::max (m_Nx - 1, 1));
    m_hy = (m_extrema->y - m_origin->y) / double (std::max (m_Ny - 1, 1));
    m_hz = (m_extrema->z - m_origin->z) / double (std::max (m_Nz - 1, 1));

    double eps = 1e-8;

    m_hx = (m_hx < eps ? 1.: m_hx);
    m_hy = (m_hy < eps ? 1.: m_hy);
    m_hz = (m_hz < eps ? 1.: m_hz);

    m_points = {};

    // Définition des points (vecteurs) de R^3 dans les trois directions
    Point p = {m_hx, m_hy, m_hz};

    for (int k = 0; k < m_Nz; ++k)
        for (int j = 0; j < m_Ny; ++j)
            for (int i = 0; i < m_Nx; ++i)
                AddPointOnDomain (*m_origin + (Point (i, j, k) * p));

    for (int k = 0; k < m_Nz; ++k)
        for (int j = 0; j < m_Ny; ++j)
            for (int i = 0; i < m_Nx; ++i)
            {
                Point* tp = GetPoint (i, j, k);

                if (m_Nx > 1)
                {
                    tp->AddPointNeighbour (GetPoint (i-1, j, k));

                    if (m_Nx > 2)
                        tp->AddPointNeighbour (GetPoint (i+1, j, k));

                }

                if (m_Ny > 1)
                {
                    tp->AddPointNeighbour (GetPoint (i, j-1, k));

                    if (m_Ny > 2)
                        tp->AddPointNeighbour (GetPoint (i, j+1, k));

                }

                if (m_Nz > 1)
                {
                    tp->AddPointNeighbour (GetPoint (i, j, k-1));

                    if (m_Nz > 2)
                        tp->AddPointNeighbour (GetPoint (i, j, k+1));

                }

            }

    for (int k = 0; k < std::max (m_Nz - 1, 1); ++k)
        for (int j = 0; j < std::max (m_Ny - 1, 1); ++j)
            for (int i = 0; i < std::max (m_Nx - 1, 1); ++i)
            {
                std::vector <Point *> list_points = {};

                // p_{i, j, k}
                list_points.push_back (GetPoint (i, j, k));

                // p_{i+1, j, k}
                list_points.push_back (GetPoint (i + 1, j, k));

                if (m_Ny > 1)
                {
                    // p_{i+1, j+1, k}
                    list_points.push_back (GetPoint (i + 1, j + 1, k));
                    // p_{i, j+1, k}
                    list_points.push_back (GetPoint (i, j + 1, k));
                }

                if (m_Nz > 1)
                {
                    // p_{i, j, k + 1}
                    list_points.push_back (GetPoint (i, j, k + 1));

                    // p_{i+1, j, k+1}
                    list_points.push_back (GetPoint (i + 1, j, k + 1));

                    // p_{i+1, j+1, k+1}
                    list_points.push_back (GetPoint (i + 1, j + 1, k + 1));

                    // p_{i, j+1, k+1}
                    list_points.push_back (GetPoint (i, j + 1, k + 1));
                }

                MakeACellFromListPoints (list_points);
            }
    return this;
}

Point* Mesh::GetPoint (int index) const
{
    return m_points.at (size_t (index));
}

Point* Mesh::GetPoint (int i, int j, int k) const
{
    return GetPoint (IndexPoints (i, j, k));
}

Cell* Mesh::GetCell (int index) const
{
    return m_cells.at (size_t (index));
}

Cell* Mesh::GetCell (int i, int j, int k) const
{
    return GetCell (IndexCells (i, j, k));
}

DIM Mesh::GetDimension () const
{
    if (m_Nz >= 1)
        return DIM_3D;
    else if (m_Ny >= 2)
        return DIM_2D;
    return DIM_1D;
}

std::vector<Point *> Mesh::GetBounds () const
{
    return {m_origin, m_extrema};
}

int Mesh::Get_Nx () const
{
    return m_Nx;
}

int Mesh::Get_Ny () const
{
    return m_Ny;
}

int Mesh::Get_Nz () const
{
    return m_Nz;
}

double Mesh::Get_hx () const
{
    return m_hx;
}

double Mesh::Get_hy () const
{
    return m_hy;
}

double Mesh::Get_hz () const
{
    return m_hz;
}

int Mesh::GetNumberOfTotalPoints () const
{
    return int (m_points.size ());
}

int Mesh::GetNumberOfCartesianPoints () const
{
    return m_Nx * m_Ny * m_Nz;
}

int Mesh::GetNumberOfTotalCells () const
{
    return int (m_cells.size ());
}

void Mesh::Print ()
{
    int G = GetNumberOfCartesianPoints ();
    int T = int (m_points.size ());
    int B = int (GetListOfIndexPoints().size());

    std::cout << "-- MESH CLASS INTERNAL PRINT ---" << std::endl;
    std::cout << "Origin                = " << *m_origin << std::endl;
    std::cout << "Extrema               = " << *m_extrema << std::endl;
    std::cout << "{N_x, Ny, Nz}         = {" << m_Nx << ", " << m_Ny << ", " << m_Nz << "}" << std::endl;
    std::cout << "{h_x, hy, hz}         = {" << m_hx << ", " << m_hy << ", " << m_hz << "}" << std::endl;
    std::cout << "Grid Points           = " << G << std::endl;
    std::cout << "Border Points         = " << B << std::endl;
    std::cout << "Grid Points on Border = " << G-T+B << std::endl << std::endl;

    std::cout << "TOTAL Points      = " << m_points.size () << std::endl;
    std::cout << "TOTAL Cells       = " << m_cells.size () << "\n" << std::endl;


    //    std::cout << "List Points : " << std::endl;

    //    for (int k = 0; k < m_Nz; ++k)
    //        for (int j = 0; j < m_Ny; ++j)
    //            for (int i = 0; i < m_Nx; ++i)
    //            {
    //                int index = IndexPoints (i, j, k);
    //                std::cout << "p_" << index << ":\t" << m_points [size_t (index)] << std::endl;
    //            }
    return;
}

Mesh * Mesh::TagPoint (int i, int j, int k, POINT_LOCATION tag)
{
    return TagPoint (IndexPoints (i, j, k), tag);
}

Mesh * Mesh::TagPoint (int index, POINT_LOCATION tag)
{
    m_points.at (size_t (index))->SetLocate (tag);
    return this;
}

Mesh* Mesh::AddPointOnBorder (Point a)
{
    auto p = new Point (a);
    p->SetLocate (ON_BORDER_OMEGA);
    p->SetGlobalIndex (int (m_points.size ()));
    m_points.push_back (p);

    return this;
}

Mesh* Mesh::AddPointOnBorder (Point* a)
{
    a->SetLocate (ON_BORDER_OMEGA);
    a->SetGlobalIndex (int (m_points.size ()));
    m_points.push_back (a);

    return this;
}

Mesh* Mesh::AddPointOnDomain (Point a, POINT_LOCATION tag)
{
    auto p = new Point (a);
    p->SetLocate (tag);
    p->SetGlobalIndex (int (m_points.size ()));
    m_points.push_back (p);
    return this;
}

Mesh* Mesh::AddPointOnDomain (Point * a, POINT_LOCATION tag)
{
    a->SetLocate (tag);
    a->SetGlobalIndex (int (m_points.size ()));
    m_points.push_back (a);
    return this;
}

Mesh* Mesh::MakeACellFromListPoints (std::vector<Point *> list)
{
    Cell* c = new Cell ();
    for (Point* p : list)
        c->AddPoint (p);

    m_cells.push_back (c);

    return this;
}

std::vector <int> Mesh::GetListOfIndexPoints (POINT_LOCATION tag)
{
    std::vector <int> indexes = {};
    for (unsigned int i = 0; i < m_points.size (); ++i)
    {
        if (m_points.at (i)->GetLocate () == tag)
            indexes.push_back (int (i));
    }
    return indexes;
}

std::vector <int> Mesh::GetListOfIndexCells (CELL_LOCATION tag)
{
    std::vector <int> indexes = {};
    for (unsigned int i = 0; i < m_cells.size (); ++i)
    {
        if (m_cells.at (i)->GetLocate ()== tag)
            indexes.push_back (int (i));
    }
    return indexes;
}

int Mesh::GetNumberOfInfosCells () const
{
    int sum = 0;
    for (size_t i = 0; i < m_cells.size (); ++i)
        sum += m_cells.at (i)->GetNumberOfInfos ();
    return sum;
}

int Mesh::IndexPoints (int i, int j, int k) const
{

    while (i >= m_Nx || i < 0)
        i += (i < 0 ? 1 : -1) * m_Nx;

    while (j >= m_Ny || j < 0)
        j += (j < 0 ? 1 : -1) * m_Ny;

    while (k >= m_Nz || k < 0)
        k += (k < 0 ? 1 : -1) * m_Nz;

    int index = (k * m_Ny + j) * m_Nx + i;

    //    printf("Call index point : %i\n", index);

    if (index < this->GetNumberOfTotalPoints () && index >= 0)
        return index;

    return 0;
}

Mesh* Mesh::RemoveThisPoint (int index)
{
    m_points.erase (m_points.begin () + index);

    return this;
}

Mesh* Mesh::RemoveThisPoint (int i, int j, int k)
{
    return RemoveThisPoint (IndexPoints (i, j, k));
}

Mesh* Mesh::RemoveThesePoints (std::vector<int> listindex)
{

    for (int index : listindex)
        RemoveThisPoint (index);

    return this;
}

Mesh* Mesh::RemoveAllNotCartesianPoints ()
{
    for (Point * p : m_points)
    {
        if (p->GetGlobalIndex () >= GetNumberOfCartesianPoints ())
            RemoveThisPoint (p->GetGlobalIndex ());
    }

    return this;
}

int Mesh::IndexCells (int i, int j, int k) const
{
    int index = (k * (m_Ny - 1) + j) * (m_Nx - 1) + i;

    if (index < this->GetNumberOfTotalCells () && index >= 0)
        return index;
    return 0;
}










