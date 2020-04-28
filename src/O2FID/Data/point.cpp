/** @file point.cpp */

#include "point.h"
#include "cell.h"

Point::Point () :
    x (0.), y (0.), z (0.),
    m_cells ({}),
    m_listNeighbours ({}),
    m_globalIndex (-1), // Il n'a pas encore d'indice global
    m_locate (ON_DOMAIN_EXTERN_OMEGA) // Par défaut le point est considéré à l'extèrieur
{}

Point::Point (const Point &p) :
    x (p.x), y (p.y), z (p.z),
    m_cells (p.m_cells),
    m_listNeighbours (p.m_listNeighbours),
    m_globalIndex (p.m_globalIndex),
    m_locate (p.m_locate)
{}

Point::Point (double a, double b, double c) :
    x (a), y (b), z (c),
    m_cells ({}),
    m_listNeighbours ({}),
    m_globalIndex (-1), // Il n'a pas encore d'indice global
    m_locate (ON_DOMAIN_EXTERN_OMEGA) // Par défaut le point est considéré à l'extèrieur
{}

Point& Point::operator= (const Point& p)
{
    x = p.x;
    y = p.y;
    z = p.z;
    m_cells  = p.m_cells;
    m_listNeighbours  = p.m_listNeighbours;
    m_globalIndex = p.m_globalIndex;
    m_locate = p.m_locate;

    return *this;
}

Point& Point::operator= (std::initializer_list<double> ilist)
{
    auto iter = ilist.begin ();
    x = (ilist.size () >= 1 ? *iter : 0.);
    iter++;
    y = (ilist.size () >= 2 ? *iter : 0.);
    iter++;
    z = (ilist.size () >= 3 ? *iter : 0.);

    return *this;
}

Point& Point::operator= (std::initializer_list<int> ilist)
{
    auto iter = ilist.begin ();
    x = double (ilist.size () >= 1 ? *iter : 0.);
    iter++;
    y = double (ilist.size () >= 2 ? *iter : 0.);
    iter++;
    z = double (ilist.size () >= 3 ? *iter : 0.);

    return *this;
}

Point::~Point ()
{}


Point* Point::SetGlobalIndex (int index)
{
    m_globalIndex = index;
    return this;
}

int Point::GetGlobalIndex () const
{
    return m_globalIndex;
}

Point* Point::SetLocate (POINT_LOCATION loc)
{
    m_locate = loc;
    return this;
}

POINT_LOCATION Point::GetLocate () const
{
    return m_locate;
}

void Point::AddPointNeighbour (Point* p)
{
    m_listNeighbours.push_back (p);
    return;
}

void Point::RemoveThisNeighbourPoint (Point* p)
{
    auto it = m_listNeighbours.begin ();

    while (it != m_listNeighbours.end ())
    {
        if (p == *it)
            it = m_listNeighbours.erase (it);
        else
            it++;
    }

    return;
}

void Point::ClearListNeighbours ()
{
    m_listNeighbours.clear ();
    return;
}

std::vector<Point*> Point::GetListNeighbours ()
{
    return m_listNeighbours;
}

void Point::SetListNeighbours (std::vector <Point *>& list)
{
    for (Point * p : list)
        if (p == nullptr)
            return;

    m_listNeighbours = list;

    return;
}

void Point::LinkToCell (Cell *c)
{
    for (Cell* tc : m_cells)
    {
        if (tc == c)
            return;
    }

    m_cells.push_back (c);

    return;
}

void Point::UnlinkToCell (Cell *c)
{
    auto it = m_cells.begin ();
    while (it != m_cells.end ())
    {
        if (*it == c)
            m_cells.erase (it);
        else
            it++;
    }

    return;
}

std::vector <Cell*> Point::GetLinkedCell () const
{
    return m_cells;
}

Point* Point::DetachFromAll ()
{
    for (Cell* c : m_cells)
        c->RemovePoint (this);

    for (Point* p : m_listNeighbours)
        p->RemoveThisNeighbourPoint (this);
    if (m_listNeighbours.size () == 2)
    {
        m_listNeighbours.at (0)->AddPointNeighbour (m_listNeighbours.at (1));
        m_listNeighbours.at (1)->AddPointNeighbour (m_listNeighbours.at (0));
    }

    m_listNeighbours.clear ();
    m_cells.clear ();

    return this;
}

std::ostream & operator<< (std::ostream &out, const Point &p)
{
    out << SPACE << p.x << " "
        << SPACE << p.y << " "
        << SPACE << p.z << " ";
    return out;
}

double EuclidianDist (const Point& a, const Point& b)
{
    double x = a.x - b.x;
    double y = a.y - b.y;
    double z = a.z - b.z;

    return std::sqrt (x * x + y * y + z * z);
}

Point operator* (double value, const Point& p)
{
    return {value * p.x, value * p.y, value * p.z};
}

Point operator* (const Point& p, double value)
{
    return value * p;
}

Point operator* (const Point& a, const Point& b)
{
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}

Point operator+ (const Point& a, const Point& b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Point operator- (const Point& a, const Point& b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

double operator| (const Point& a, const Point& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Point operator+ (double value, const Point& p)
{
    return {value + p.x, value + p.y, value + p.z};
}

Point operator+ (const Point& p, double value)
{
    return value + p;
}

Point operator- (double value, const Point& p)
{
    return {value - p.x, value - p.y, value - p.z};
}

Point operator- (const Point& p, double value)
{
    return (-1 * value) + p;
}

Point operator/ (double value, const Point& p)
{
    return {value / p.x, value / p.y, value / p.z};
}

Point operator/ (const Point& p, double value)
{
    return (1. / value) * p;
}

bool operator< (const Point& a, const Point& b)
{
    return ((a.x < b.x) && (a.y < b.y) && (a.z < b.z));
}

bool operator> (const Point& a, const Point& b)
{
    return ((a.x > b.x) && (a.y > b.y) && (a.z > b.z));
}

bool operator<= (const Point& a, const Point& b)
{
    return ((a.x <= b.x) && (a.y <= b.y) && (a.z <= b.z));
}

bool operator>= (const Point& a, const Point& b)
{
    return ((a.x >= b.x) && (a.y >= b.y) && (a.z >= b.z));
}

bool operator== (const Point& a, const Point& b)
{
    double eps = 1e-10;
    bool bx = (fabs(a.x - b.x) < eps);
    bool by = (fabs(a.y - b.y) < eps);
    bool bz = (fabs(a.z - b.z) < eps);

    return (bx && by && bz);
}

bool operator!= (const Point& a, const Point& b)
{
    double eps = 1e-10;
    bool bx = (fabs(a.x - b.x) > eps);
    bool by = (fabs(a.y - b.y) > eps);
    bool bz = (fabs(a.z - b.z) > eps);

    return (bx && by && bz);
}
