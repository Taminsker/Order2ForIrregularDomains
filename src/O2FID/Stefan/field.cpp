/** @file field.cpp */

#include "field.h"

Field::Field() :
    W(std::vector<Point*>()),
    Normals(std::vector<Point*>()),
    GradTemperature(std::vector<Point*>()),
    GradPhi(std::vector<Point*>())
{}

Field::Field(const Field& f)
{
    W.resize (f.W.size ());
    for (size_t i = 0; i < f.W.size (); ++i)
    {
        W.at (i) = new Point (*f.W.at (i));
    }

    Normals.resize (f.Normals.size ());
    for (size_t i = 0; i < f.Normals.size (); ++i)
    {
       Normals.at (i) = new Point (*f.Normals.at (i));
    }

    GradPhi.resize (f.GradPhi.size ());
    for (size_t i = 0; i < f.GradPhi.size (); ++i)
    {
       GradPhi.at (i) = new Point (*f.GradPhi.at (i));
    }

    GradTemperature.resize (f.GradTemperature.size ());
    for (size_t i = 0; i < f.GradPhi.size (); ++i)
    {
       GradTemperature.at (i) = new Point (*f.GradTemperature.at (i));
    }
}


Field::~Field ()
{
    AutoClearVector(&W);
    AutoClearVector(&Normals);
    AutoClearVector(&GradTemperature);
    AutoClearVector(&GradPhi);
}

Field& Field::operator= (const Field& f)
{
    W.resize (f.W.size ());
    for (size_t i = 0; i < f.W.size (); ++i)
    {
        W.at (i) = new Point (*f.W.at (i));
    }

    Normals.resize (f.Normals.size ());
    for (size_t i = 0; i < f.Normals.size (); ++i)
    {
       Normals.at (i) = new Point (*f.Normals.at (i));
    }

    GradPhi.resize (f.GradPhi.size ());
    for (size_t i = 0; i < f.GradPhi.size (); ++i)
    {
       GradPhi.at (i) = new Point (*f.GradPhi.at (i));
    }

    GradTemperature.resize (f.GradTemperature.size ());
    for (size_t i = 0; i < f.GradPhi.size (); ++i)
    {
       GradTemperature.at (i) = new Point (*f.GradTemperature.at (i));
    }

    return *this;
}
