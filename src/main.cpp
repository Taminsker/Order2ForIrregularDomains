#include <iostream>
#include "Chaleur1D.h"

int main(int argc, char* argv[])
{
    auto v = mesh (0, 1, 10);
    auto m = matrix (0, 1, 10, 1, 1);

    std::cout << v << std::endl;
    return 0;
}
