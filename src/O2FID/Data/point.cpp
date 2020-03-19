/** @file point.cpp */

#include "point.h"

Point::Point () :
    x (0.), y (0.), z (0.),
    m_locate (ON_DOMAIN_EXTERN_OMEGA), // Par défaut le point est considéré à l'extèrieur
    m_globalIndex (-1) // Il n'a pas encore d'indice global
{}

Point::Point (const Point &p) :
    x (p.x), y (p.y), z (p.z),
    m_cells (p.m_cells),
    m_locate (p.m_locate),
    m_globalIndex (p.m_globalIndex),
    m_listNeighbours (p.m_listNeighbours)
{}

Point::Point (double a, double b, double c) :
    x (a), y (b), z (c),
    m_locate (ON_DOMAIN_EXTERN_OMEGA), // Par défaut le point est considéré à l'extèrieur
    m_globalIndex (-1) // Il n'a pas encore d'indice global
{}

Point::~Point () {}

void Point::SetGlobalIndex (int index)
{
    m_globalIndex = index;
    return;
}

int Point::GetGlobalIndex ()
{
    return m_globalIndex;
}

void Point::SetLocate (POINT_LOCATION loc)
{
    m_locate = loc;
    return;
}

POINT_LOCATION Point::GetLocate () const
{
    return m_locate;
}

void Point::AddPointNeighbour (Point p)
{
    m_listNeighbours.push_back (p);
    return;
}

void Point::ClearListNeighboors ()
{
    m_listNeighbours.clear ();
    return;
}

std::vector <Point> Point::GetListNeighbours ()
{
    return m_listNeighbours;
}

bool Point::operator== (const Point &p)
{
    double eps = 1e-10;
    if (    fabs (x - p.x) < eps &&
            fabs (y - p.y) < eps &&
            fabs (z - p.z) < eps)
        return true;
    return false;
}

Point Point::operator+ (const Point &a)
{
    return Point (x + a.x, y + a.y, z + a.z);
}

Point Point::operator- (const Point &a)
{
    return Point (x - a.x, y - a.y, z - a.z);
}

Point Point::operator* (double v)
{
    return Point (v * x, v * y, v * z);
}

Point Point::operator/ (double v)
{
    if (std::abs(v) < 1e-10)
        return *this;
    return Point (x / v, y / v, z / v);
}

Point Point::operator * (int v)
{
    return Point (v * x, v * y, v * z);
}
std::ostream & operator<< (std::ostream &out, const Point &p)
{
    out << SPACE << p.x << " "
        << SPACE << p.y << " "
        << SPACE << p.z << " ";
    return out;
}
