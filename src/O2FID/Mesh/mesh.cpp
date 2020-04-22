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
    std::cout << "# Mesh builder" << std::endl;

    if (m_Nx == 0 || m_Ny == 0 || m_Nz == 0) // Pas de points
    {
        std::cout << INDENT << "ERROR Mesh:: triplet for {Nx, Ny, Nz} "
                  << m_Nx << " " << m_Ny << " " << m_Nz << std::endl;
        exit(0);
    }

    if (*m_origin == *m_extrema) // /!\ domaine réduit à un point
    {
        std::cout << INDENT << "ERROR Mesh:: vector bounds for "
                  << "[" << m_origin->x << ", " << m_extrema->x << "] x "
                  << "[" << m_origin->y << ", " << m_extrema->y << "] x "
                  << "[" << m_origin->z << ", " << m_extrema->z << "]"
                  << std::endl;
        m_Nx = m_Ny = m_Nz = 0;
        exit(0);
    }

    int request = m_Nx * m_Ny * m_Nz;
    std::cout << INDENT << "You request " << request << " points." << std::endl;

    if (request > 1e6)
    {
        char answer = 't';
        while (answer != 'y' && answer != 'n')
        {
            std::cout << "\r" << INDENT << "Are you sure to continue ? (y/n) : " << std::flush;
            std::cin >> answer;
        }

        std::cout << std::endl;
        if (answer == 'n')
        {
            std::cout << "ABORT" << std::endl;
            exit(0);
        }
    }

    m_hx = (m_extrema->x - m_origin->x) / double (std::max (m_Nx - 1, 1));
    m_hy = (m_extrema->y - m_origin->y) / double (std::max (m_Ny - 1, 1));
    m_hz = (m_extrema->z - m_origin->z) / double (std::max (m_Nz - 1, 1));

    double eps = 1e-20;

    m_hx = (m_hx < eps ? 0.: m_hx);
    m_hy = (m_hy < eps ? 0.: m_hy);
    m_hz = (m_hz < eps ? 0.: m_hz);

    bool b1 = (m_hx < eps ? true : false);
    bool b2 = (m_hy < eps ? true : false);
    bool b3 = (m_hz < eps ? true : false);

    if (b1 && b2 && b3)
    {
        std::cout << INDENT << "ERROR Mesh:: triplet for {hx, hy, hz} "
                  << m_hx << " " << m_hy << " " << m_hz << std::endl;
        exit(1);
    }

    m_points = {};

    // Définition des points (vecteurs) de R^3 dans les trois directions
    Point p = {m_hx, m_hy, m_hz};

    for (int k = 0; k < m_Nz; ++k)
        for (int j = 0; j < m_Ny; ++j)
            for (int i = 0; i < m_Nx; ++i)
                AddPointOnDomain (*m_origin + (Point (i, j, k) * p));

    std::cout << INDENT << "Cartesian grid points have been added." << std::endl;

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

    std::cout << INDENT << "The links between the grid points have been created." << std::endl;

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

    std::cout << INDENT << "The cells were created." << std::endl;

    std::cout << std::endl;

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

Point* Mesh::GetSnakePoint (int i, int j, int k) const
{

    if (Remainder (k, 2) == 0) // k pair
    {
        //        Parcours :
        //        16      17      18      19    --> j = 4
        //        15      14      13      12    --> j = 3
        //        8       9       10      11    --> j = 2
        //        7       6       5       4     --> j = 1
        //        0       1       2       3     --> j = 0
        //        i = 0   i = 1   i = 2   i = 3
        //        i croissant sur j pair, i décroissant sur j impair, j croissant sur tout i

        if (Remainder (j, 2) == 0)
            return GetPoint (i, j, k);

        return GetPoint (m_Nx - i-1, j, k);
    }
    else
    {
        //        Parcours :
        //        3       2       1        0    --> j = 4
        //        4       5       6        7    --> j = 3
        //        11      10      9        8    --> j = 2
        //        12      13      14       15   --> j = 1
        //        19      18      17       16   --> j = 0
        //        i = 0   i = 1   i = 2    i = 3
        //        j décroissant sur tout i, i décroissant sur Ny-j pair et croissant sur Ny-j impair

        if (Remainder (m_Ny - j, 2) == 0)
            return GetPoint (m_Nx - 1 - i, m_Ny - 1 - j, k);
        return GetPoint (i, m_Ny - 1 - j, k);
    }
}

Point* Mesh::GetSnakePoint (int idx) const
{
    int i = Remainder (idx, m_Nx);
    int j = Remainder ((idx - i) / m_Nx, m_Ny);
    int k = Quotient ((idx - i) / m_Nx - j, m_Ny);

    return GetSnakePoint (i, j, k);
}

LocalIndexes Mesh::GetLocalIndexesOfPoint (int idx) const
{
    LocalIndexes r;
    r.i = Remainder (idx, m_Nx);
    r.j = Remainder ((idx - r.i) / m_Nx, m_Ny);
    r.k = Quotient ((idx - r.i) / m_Nx - r.j, m_Ny);

    return r;
}

int Mesh::GetGlobalIndexOfPoint (int i, int j, int k) const
{
    i = Remainder (i, m_Nx);
    j = Remainder (j, m_Ny);
    k = Remainder (k, m_Nz);

    return (k * m_Ny + j) * m_Nx + i;
}

Cell* Mesh::GetCell (int index) const
{
    return m_cells.at (size_t (index));
}

Cell* Mesh::GetCell (int i, int j, int k) const
{
    return GetCell (IndexCells (i, j, k));
}

LocalIndexes Mesh::GetLocalIndexesOfCell (int idx) const
{
    LocalIndexes r;
    r.i = Remainder (idx, m_Nx - 1);
    r.j = Remainder ((idx - r.i) / (m_Nx - 1), m_Ny - 1);
    r.k = Quotient ((idx - r.i) / (m_Nx - 1) - r.j, m_Ny - 1);

    return r;
}

int Mesh::GetGlobalIndexOfCell (int i, int j, int k) const
{
    i = Remainder (i, m_Nx - 1);
    j = Remainder (j, m_Ny - 1);
    k = Remainder (k, m_Nz - 1);

    return (k * (m_Ny - 1) + j) * (m_Nx - 1) + i;
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

Mesh* Mesh::Print ()
{
    int G = GetNumberOfCartesianPoints ();
    int T = int (m_points.size ());
    int B = int (GetListOfIndexPoints().size());

    std::cout << "# Mesh print" << std::endl;
    std::cout << INDENT << "Origin                = " << *m_origin << std::endl;
    std::cout << INDENT << "Extrema               = " << *m_extrema << std::endl;
    std::cout << INDENT << "{N_x, Ny, Nz}         = {" << m_Nx << ", " << m_Ny << ", " << m_Nz << "}" << std::endl;
    std::cout << INDENT << "{h_x, hy, hz}         = {" << m_hx << ", " << m_hy << ", " << m_hz << "}" << std::endl;
    std::cout << INDENT << "Grid Points           = " << G << std::endl;
    std::cout << INDENT << "Border Points         = " << B << std::endl;
    std::cout << INDENT << "Grid Points on Border = " << G-T+B << std::endl << std::endl;

    std::cout << INDENT << "TOTAL Points      = " << m_points.size () << std::endl;
    std::cout << INDENT << "TOTAL Cells       = " << m_cells.size () << "\n" << std::endl;


    //    std::cout << "List Points : " << std::endl;

    //    for (int k = 0; k < m_Nz; ++k)
    //        for (int j = 0; j < m_Ny; ++j)
    //            for (int i = 0; i < m_Nx; ++i)
    //            {
    //                int index = IndexPoints (i, j, k);
    //                std::cout << "p_" << index << ":\t" << m_points [size_t (index)] << std::endl;
    //            }

    std::cout << std::endl;

    return this;
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
    i = Remainder (i, m_Nx);
    j = Remainder (j, m_Ny);
    k = Remainder (k, m_Nz);

    int index = (k * m_Ny + j) * m_Nx + i;

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

    i = Remainder (i, m_Nx - 1);
    j = Remainder (j, m_Ny - 1);
    k = Remainder (k, m_Nz - 1);

    int index = (k * (m_Ny - 1) + j) * (m_Nx - 1) + i;

    if (index < this->GetNumberOfTotalCells () && index >= 0)
        return index;
    return 0;
}










