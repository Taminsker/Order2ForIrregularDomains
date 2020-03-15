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
        return;
    }

    if (m_origin == m_extrema)
    {
        std::cout << "ERROR Mesh:: vector bounds for "
                  << "[" << m_origin.x << ", " << m_extrema.x << "] x ["
                  << "[" << m_origin.y << ", " << m_extrema.y << "] x ["
                  << "[" << m_origin.z << ", " << m_extrema.z << "]"
                  << std::endl;
        m_Nx = m_Ny = m_Nz = 0;
        return;
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

    m_points = std::vector<Point> (static_cast<unsigned int> (m_Nx * m_Ny * m_Nz));

    for (int i = 0; i < m_Nx; ++i)
    {
        double x = m_origin.x + i * m_hx;

        for (int j = 0; j < m_Ny; ++j)
        {
            double y = m_origin.y + j * m_hy;

            for (int k = 0; k < m_Nz; ++k)
            {
                double z = m_origin.z + k * m_hz;

                m_points [static_cast<unsigned int> (Index (i, j, k))] = Point (x, y, z);
            }
        }
    }

    return;
}

Point Mesh::operator() (int index) const
{
    return m_points [static_cast<unsigned int>(index)];
}

Point Mesh::operator()(int i, int j, int k) const
{
    return m_points [static_cast<unsigned int>(Index (i, j, k))];
}

int Mesh::GetDimension ()
{
    return m_dim;
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
                int index = Index (i, j, k);
                std::cout << "p_" << index << ":\t" << m_points [static_cast<unsigned int>(index)] << std::endl;
            }
    return;
}

int Mesh::AddPointOnBorder (Point a)
{
    a.SetLocate (ON_BORDER_OMEGA);
    m_points.push_back (a);
    return m_points.size () - 1;
}

int Mesh::AddPointOnDomain (Point a)
{
    a.SetLocate (ON_DOMAIN_OMEGA);
    m_points.push_back (a);
    return m_points.size () - 1;
}

std::vector <int> Mesh::GetListOfIndexPoints ()
{
    std::vector <int> indexes = {};
    for (unsigned int i = 0; i < m_points.size (); ++i)
    {
        if (m_points.at (i).GetLocate () == ON_BORDER_OMEGA)
            indexes.push_back (int (i));
    }
    return indexes;
}

int Mesh::Index (int i, int j, int k) const
{
    return int ((k * m_Ny + j) * m_Nx + i);
}
