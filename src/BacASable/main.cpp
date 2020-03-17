#include "../O2FID/O2FID.h"

#include <iostream>

double phi (Point a, double t = 0); // fonction levelset
double f_laplace (Point a, double t = 0.); // fonction de second membre
double g_laplace (Point a, double t = 0.); // fonction des conditions aux limites

double f_chaleur (Point a, double t); // fonction de second membre equation chaleur
double g_chaleur (Point a, double t); // fonction des conditions aux limites equation chaleur

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    Mesh * mesh = new Mesh ();
    mesh->Set_Nx (100);
    mesh->Set_Ny (100);
    mesh->Set_Nz (100);
    mesh->SetBounds (Point (1, 1, 1), Point (1, 1, 1));
    mesh->Build ();
    mesh->Print ();

    Vector u (mesh->GetNumberOfTotalPoints ());
    u.setOnes ();

    Writer writer (mesh);
    writer.SetVectorAnalytical (&u);

    writer.SetFilename ("Test");
//    writer.WriteNow ();

    return 0;
}

double phi (Point a, double t)
{
    (void)t;
    return 0. * (a.x + a.y + a.z);
}

double f_laplace (Point a, double t)
{
    (void)t;
    return 0. * (a.x + a.y + a.z);
}

double g_laplace (Point a, double t)
{
    (void)t;
    return 0. * (a.x + a.y + a.z);
}

double f_chaleur (Point a, double t)
{
    return 0. * (a.x + a.y + a.z) + t;
}

double g_chaleur (Point a, double t)
{
    return 0. * (a.x + a.y + a.z) + t;
}
