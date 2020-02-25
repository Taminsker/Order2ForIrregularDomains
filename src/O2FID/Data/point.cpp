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

bool Point::operator== (const Point &p)
{
    double eps = 1e-6;
    if (    fabs (x - p.x) < eps &&
            fabs (y - p.y) < eps &&
            fabs (z - p.z) < eps)
        return true;
    return false;
}

std::ostream & operator<< (std::ostream &out, const Point &p)
{
    out << SPACE << p.x << " "
        << SPACE << p.y << " "
        << SPACE << p.z << " ";
    return out;
}
