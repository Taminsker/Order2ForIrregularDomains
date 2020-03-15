/** @file point.cpp */

#include "point.h"

Point::Point () :
    x (0.), y (0.), z (0.)
{}

Point::Point (const Point &p) :
    x (p.x), y (p.y), z (p.z)
{}

Point::Point (double a, double b, double c) :
    x (a), y (b), z (c)
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

void Point::SetLocate (int loc)
{
    m_locate = loc;
    return;
}

int Point::GetLocate () const
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
    double eps = 1e-6;
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

std::ostream & operator<< (std::ostream &out, const Point &p)
{
    out << SPACE << p.x << " "
        << SPACE << p.y << " "
        << SPACE << p.z << " ";
    return out;
}
