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

    double time_final = 0.248;
    auto Nx = {81, 161, 321};
    auto Ny = {81, 161, 321};
//    auto Nx = {321};
//    auto Ny = {321};
    auto No = {0};

    Make (&min, &max, time_final, BackwardEuler1, Nx, Ny, No);
    Make (&min, &max, time_final, CrankNicolson, Nx, Ny, No);
    Make (&min, &max, time_final, BackwardEuler2, Nx, Ny, No);

    return 0;
}
