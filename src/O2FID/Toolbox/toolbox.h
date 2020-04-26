#ifndef TOOLBOX_H
#define TOOLBOX_H

#define INDENT "-->\t"

#include <iostream>
#include <vector>
#include <math.h>

#include "../Data/datatypedefinitions.h"

/*!
 *  \addtogroup Outils
 *  @{
 */

int Remainder (int dividend, int divisor);
int Quotient (int dividend, int divisor);

std::vector<double> Order (std::vector<double> err, std::vector<double> h);

template<typename T>
std::ostream & operator<< (std::ostream &out, const std::vector<T> vec)
{
    for (size_t i = 0; i < vec.size (); ++i)
        out << vec.at (i) << " " << std::flush;
    return out;
}

class Point;
template <typename T>
std::vector<T*>* AutoClearVector(std::vector<T*>* vec)
{
    for (size_t i = 0; i < vec->size (); ++i)
    {
        delete vec->at (i);
        vec->at (i) = nullptr;
    }

    vec->clear ();

    return vec;
}

template <typename T>
std::vector<T> operator* (T value, std::vector<T> vec)
{
    std::vector<T> R(vec.size (), T(0));

    for (size_t i = 0; i < vec.size (); ++i)
    {
        R.at (i) = value * vec.at (i);
    }

    return R;
}

class Mesh;

void Extrapole (Mesh* mesh, Vector* vec);

void Extrapole (Mesh* mesh, std::vector<Point*>* vec);


/** @} */


#endif // TOOLBOX_H
