#ifndef TOOLBOX_H
#define TOOLBOX_H

#define INDENT "-->\t"

#include <iostream>
#include <vector>
#include <math.h>

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

/** @} */


#endif // TOOLBOX_H
