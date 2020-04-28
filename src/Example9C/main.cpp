#include "../O2FID/O2FID.h"

#include "headers.h"

#include <iostream>
#include <cmath>
#include "Eigen/SparseCore"

#define SPACER std::left << std::setw(20)




int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "            EXAMPLE 9 - O2FID            " << std::endl;
    std::cout << "-----------------------------------------" << std::endl;


    Point min = {0., 0., 0.};
    Point max = {0.5, 0.5, 0.5};

    int Nt = 11;
    auto Nx = {26, 51, 101};

    Make (&min, &max, Nt, BackwardEuler1, Nx, Nx, Nx);
    Make (&min, &max, Nt, CrankNicolson, Nx, Nx, Nx);
    Make (&min, &max, Nt, BackwardEuler2, Nx, Nx, Nx);

    return 0;
}
