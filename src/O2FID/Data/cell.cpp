#include "cell.h"
#include "point.h"

Cell::Cell() :
    m_cell_n (nullptr), m_cell_s (nullptr),
    m_cell_o (nullptr), m_cell_e (nullptr)
{}

Cell::~Cell ()
{

}

void Cell::SetCellNorth (Cell *cell)
{
    m_cell_n = cell;
    return;
}

void Cell::SetCellSouth (Cell *cell)
{
    m_cell_s = cell;
    return;
}

void Cell::SetCellWest (Cell *cell)
{
    m_cell_o = cell;
    return;
}

void Cell::SetCellEast (Cell *cell)
{
    m_cell_e = cell;
    return;
}

void Cell::AddPoint(Point &p)
{
    m_points.push_back (&p);
    return;
}

Cell * Cell::GetCellNorth ()
{
    return m_cell_n;
}

Cell * Cell::GetCellSouth ()
{
    return m_cell_s;
}

Cell * Cell::GetCellWest ()
{
    return m_cell_o;
}

Cell * Cell::GetCellEast ()
{
    return m_cell_e;
}

Point Cell::GetBarycenter ()
{
    Point p = Point ();
    int N = int (m_points.size ());
    for (int i = 0; i < N; ++i)
    {
        auto a = m_points.at (static_cast<unsigned int> (i));
        p.x += a->x;
        p.y += a->y;
        p.z += a->z;
    }

    return p / double (N);
}
