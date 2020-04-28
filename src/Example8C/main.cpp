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
    std::cout << "            EXAMPLE 8 - O2FID            " << std::endl;
    std::cout << "-----------------------------------------" << std::endl;


    Point min = {-1., -1, 0};
    Point max = {1., 1, 0};

    int Nt = 21;
    auto Nx = {81, 161, 321};
    auto Ny = {81, 161, 321};
    auto No = {0, 0, 0};

    Make (&min, &max, Nt, BackwardEuler1, Nx, Ny, No);
    Make (&min, &max, Nt, CrankNicolson, Nx, Ny, No);
    Make (&min, &max, Nt, BackwardEuler2, Nx, Ny, No);

    return 0;
}
