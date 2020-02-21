#include "../O2FID/O2FID.h"

#include <iostream>

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    auto v = Mesh ();

    v.SetBounds (Point(), Point (1, 1, 0));
    v.Set_Nx(3);
    v.Set_Ny (3);
    v.Set_Nz (4);
    v.Build ();


    v.Print ();

    std::cout << "courcou" << std::endl;
    return 0;
}
