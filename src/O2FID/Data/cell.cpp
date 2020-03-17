#include "cell.h"
#include "point.h"

Cell::Cell()
{}

Cell::~Cell ()
{}

void Cell::AddPoint (Point *p)
{
    for (size_t i = 0; i < m_points.size (); ++i)
        if (m_points.at (i) == p)
            return;

    m_points.push_back (p);
    return;
}

void Cell::RemovePoint (Point *p)
{
    for (size_t i = 0; i < m_points.size (); ++i)
        if (m_points.at (i) == p)
            m_points.erase (m_points.begin () + long (i));
    return;
}

CELL_LOCATION Cell::GetLocate () const
{
    CELL_LOCATION loc = IN_DOMAIN_EXTERN_OMEGA;
    if (m_points.size () >= 1)
        loc = IN_DOMAIN_INTERN_OMEGA;
    else
        return loc;

    for (size_t i = 1; i < m_points.size (); ++i)
    {
        POINT_LOCATION loc_p = m_points.at (i)->GetLocate ();

        bool boolean1 = (loc_p == ON_DOMAIN_EXTERN_OMEGA &&
                         loc == IN_DOMAIN_INTERN_OMEGA);
        bool boolean2 = (loc_p == ON_DOMAIN_INTERN_OMEGA &&
                         loc == IN_DOMAIN_EXTERN_OMEGA);

        if (boolean1 || boolean2)
        {
            loc = IN_MIX_DOMAIN;
            break;
        }
    }
    return loc;
}

int Cell::GetType () const
{
    switch (m_points.size ()) {
    default:
    case 0:
        return 0; // VTK_EMPTY_CELL = 0;
    case 1:
        return 1; // VTK_VERTEX = 1;
    case 2:
        return 3; // VTK_LINE = 3;
    case 3:
        return 5; // VTK_TRIANGLE = 5;
    case 4:
        return 9; // VTK_QUAD = 9;
    case 5:
        return 10; // VTK_TETRA = 10;
    case 8:
        return 12; // VTK_HEXAHEDRON = 12;
    }
}

int Cell::GetNumberOfInfos () const
{
    return 1 + int (m_points.size ());
}

std::ostream& operator<< (std::ostream &out, const Cell &c)
{
    out << c.m_points.size () << " ";
    for (size_t i = 0; i < c.m_points.size (); ++i)
    {
        out << c.m_points.at (i)->GetGlobalIndex () << " ";
    }

    out << std::endl;

    return out;
}
