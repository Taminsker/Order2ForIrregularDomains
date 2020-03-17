/** @file mesh.cpp */

#include "mesh.h"

Mesh::Mesh () :
    m_hx (0.), m_hy (0.), m_hz (0.),
    m_Nx (1), m_Ny (1), m_Nz (1)
{}

Mesh::~Mesh ()
{}

void Mesh::Set_Nx (int Nx)
{
    m_Nx = std::abs (Nx);
    return;
}

void Mesh::Set_Ny (int Ny)
{
    m_Ny = std::abs (Ny);
    return;
}

void Mesh::Set_Nz (int Nz)
{
    m_Nz = std::abs (Nz);
    return;
}

void Mesh::SetBounds (Point Origin, Point Extrema)
{
    m_origin = Origin;
    m_extrema = Extrema;
    return;
}

void Mesh::Build ()
{
    if (m_Nx == 0 || m_Ny == 0 || m_Nz == 0)
    {
        std::cout << "ERROR Mesh:: triplet for {Nx, Ny, Nz} "
                  << m_Nx << " " << m_Ny << " " << m_Nz << std::endl;
        exit(0);
    }

    if (m_origin == m_extrema)
    {
        std::cout << "ERROR Mesh:: vector bounds for "
                  << "[" << m_origin.x << ", " << m_extrema.x << "] x "
                  << "[" << m_origin.y << ", " << m_extrema.y << "] x "
                  << "[" << m_origin.z << ", " << m_extrema.z << "]"
                  << std::endl;
        m_Nx = m_Ny = m_Nz = 0;
        exit(0);
    }

    m_hx = (m_extrema.x - m_origin.x) / double (std::max (m_Nx - 1, 1));
    m_hy = (m_extrema.y - m_origin.y) / double (std::max (m_Ny - 1, 1));
    m_hz = (m_extrema.z - m_origin.z) / double (std::max (m_Nz - 1, 1));

    double eps = 1e-8;
    if (m_hx < eps)
        m_Nx = 1;
    if (m_hy < eps)
        m_Ny = 1;
    if (m_hz < eps)
        m_Nz = 1;

    m_points = std::vector<Point> ();

    Point p_x = Point (m_hx);
    Point p_y = Point (0., m_hy);
    Point p_z = Point (0., 0., m_hz);

    for (int i = 0; i < m_Nx; ++i)
        for (int j = 0; j < m_Ny; ++j)
            for (int k = 0; k < m_Nz; ++k)
                AddPointOnDomain (m_origin + p_x * i + p_y * j + p_z * k);

    for (int i = 0; i < std::max (m_Nx - 1, 1); ++i)
        for (int j = 0; j < std::max (m_Ny - 1, 1); ++j)
            for (int k = 0; k < std::max (m_Nz - 1, 1); ++k)
            {
                Cell c;

                // p_{i, j, k}
                Point * p = &m_points.at (size_t (IndexPoints (i, j, k)));
                c.AddPoint (p);

                // p_{i+1, j, k}
                p = &m_points.at (size_t (IndexPoints (i + 1, j, k)));
                c.AddPoint (p);

                if (m_Ny > 1)
                {
                    // p_{i+1, j+1, k}
                    p = &m_points.at (size_t (IndexPoints (i + 1, j + 1, k)));
                    c.AddPoint (p);

                    // p_{i, j+1, k}
                    p = &m_points.at (size_t (IndexPoints (i, j + 1, k)));
                    c.AddPoint (p);
                }

                if (m_Nz > 1)
                {
                    // p_{i, j, k + 1}
                    p = &m_points.at (size_t (IndexPoints (i, j, k + 1)));
                    c.AddPoint (p);

                    // p_{i+1, j, k+1}
                    p = &m_points.at (size_t (IndexPoints (i + 1, j, k + 1)));
                    c.AddPoint (p);

                    // p_{i+1, j+1, k+1}
                    p = &m_points.at (size_t (IndexPoints (i + 1, j + 1, k + 1)));
                    c.AddPoint (p);

                    // p_{i, j+1, k+1}
                    p = &m_points.at (size_t (IndexPoints (i, j + 1, k + 1)));
                    c.AddPoint (p);
                }
                m_cells.push_back (c);
            }
    return;
}

Point Mesh::operator() (int index) const
{
    return m_points [static_cast<unsigned int>(index)];
}

Point Mesh::operator()(int i, int j, int k) const
{
    return m_points [static_cast<unsigned int>(IndexPoints (i, j, k))];
}

Point Mesh::GetPoint (int index) const
{
    return m_points [static_cast<unsigned int>(index)];
}

Point Mesh::GetPoint(int i, int j, int k) const
{
    return m_points [static_cast<unsigned int>(IndexPoints (i, j, k))];
}

Point* Mesh::GetPointPtr (int index)
{
    return &m_points [static_cast<unsigned int>(index)];
}

Point* Mesh::GetPointPtr (int i, int j, int k)
{
    return &m_points [static_cast<unsigned int>(IndexPoints (i, j, k))];
}


Cell Mesh::GetCell (int index) const
{
    return m_cells [static_cast<unsigned int>(index)];
}

Cell Mesh::GetCell (int i, int j, int k) const
{
    return m_cells [static_cast<unsigned int>(IndexCells (i, j, k))];
}

Cell* Mesh::GetCellPtr (int index)
{
    return &m_cells [static_cast<unsigned int>(index)];
}

Cell* Mesh::GetCellPtr (int i, int j, int k)
{
    return &m_cells [static_cast<unsigned int>(IndexCells (i, j, k))];
}

DIM Mesh::GetDimension ()
{
    if (m_Nz >= 1)
        return DIM_3D;
    else if (m_Ny >= 2)
        return DIM_2D;

    return DIM_1D;
}

std::vector<Point> Mesh::GetBounds () const
{
    std::vector<Point> r = std::vector<Point> (2);
    r [0] = m_origin;
    r [1] = m_extrema;

    return r;
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

int Mesh::GetNumberOfTotalCells () const
{
    return int (m_cells.size ());
}

void Mesh::Print () const
{
    std::cout << "-- MESH CLASS INTERNAL PRINT ---" << std::endl;
    std::cout << "Origin  = " << m_origin << std::endl;
    std::cout << "Extrema = " << m_extrema << std::endl;
    std::cout << "N_x     = " << m_Nx << std::endl;
    std::cout << "N_y     = " << m_Ny << std::endl;
    std::cout << "N_z     = " << m_Nz << std::endl;
    std::cout << "TOTAL   = " << m_points.size () << "\n" << std::endl;

    std::cout << "List Points : " << std::endl;

    for (int k = 0; k < m_Nz; ++k)
        for (int j = 0; j < m_Ny; ++j)
            for (int i = 0; i < m_Nx; ++i)
            {
                int index = IndexPoints (i, j, k);
                std::cout << "p_" << index << ":\t" << m_points [static_cast<unsigned int>(index)] << std::endl;
            }
    return;
}

void Mesh::TagPoint (int i, int j, int k, POINT_LOCATION tag)
{
    return TagPoint (IndexPoints (i, j, k), tag);
}

void Mesh::TagPoint (int index, POINT_LOCATION tag)
{
    std::cout << m_points.size () << std::endl;
    Point p = m_points.at (size_t (index));
    p.SetLocate (tag);
    m_points [size_t (index)] = p;

    return;
}

void Mesh::AddPointOnBorder (Point a)
{
    a.SetLocate (ON_BORDER_OMEGA);
    a.SetGlobalIndex (int (m_points.size ()));
    m_points.push_back (a);
    return;
}

void Mesh::AddPointOnDomain (Point a, POINT_LOCATION tag)
{
    a.SetLocate (tag);
    a.SetGlobalIndex (int (m_points.size ()));
    m_points.push_back (a);
    return;
}

std::vector <int> Mesh::GetListOfIndexPoints (POINT_LOCATION tag)
{
    std::vector <int> indexes = {};
    for (unsigned int i = 0; i < m_points.size (); ++i)
    {
        if (m_points.at (i).GetLocate () == tag)
            indexes.push_back (int (i));
    }
    return indexes;
}

std::vector <int> Mesh::GetListOfIndexCells (CELL_LOCATION tag)
{
    std::vector <int> indexes = {};
    for (unsigned int i = 0; i < m_cells.size (); ++i)
    {
        if (m_cells.at (i).GetLocate ()== tag)
            indexes.push_back (int (i));
    }
    return indexes;
}

int Mesh::GetNumberOfInfosCells () const
{
    int sum = 0;
    for (size_t i = 0; i < m_cells.size (); ++i)
        sum += m_cells.at (i).GetNumberOfInfos ();
    return sum;
}
int Mesh::IndexPoints (int i, int j, int k) const
{
    int index = (k * m_Ny + j) * m_Nx + i;
    if (index < this->GetNumberOfTotalPoints ())
        return index;
    return 0;
}

int Mesh::IndexCells (int i, int j, int k) const
{
    int index = (k * (m_Ny - 1) + j) * (m_Nx - 1) + i;
    if (index < this->GetNumberOfTotalCells ())
        return index;
    return 0;
}
