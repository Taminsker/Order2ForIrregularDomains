#include "../O2FID/O2FID.h"

#include <iostream>

Point phi (double theta, double t = 0); // fonction levelset
Point phi_prim (double theta, double t = 0); // fonction levelset

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
    mesh->Set_Nz (1);
    mesh->SetBounds (new Point (-1, 0), new Point (1, 3, 0));
    mesh->Build ();
    mesh->Print ();

    //    for (int i = 0; i < mesh->GetNumberOfTotalPoints (); ++i)
    //    {
    //        auto p = mesh->GetPoint (i);
    //        auto voisins = p->GetListNeighbours ();

    //        std::cout << "voisins de " << p->GetGlobalIndex () << " : " << std::flush;


    //        for (int j = 0; j < voisins.size (); ++j)
    //        {
    //            std::cout << voisins.at (j)->GetGlobalIndex ()<< " " << std::flush;
    //        }
    //        std::cout << std::endl;
    //    }




    Vector phi_value = FunToVec (mesh, phi, phi_prim);
    auto list = MakeListOfIndexPoints (mesh, &phi_value);

    mesh->Print ();
    //    for (int i : list)
    //    {
    //        //        std::cout << "nvx point : " << *mesh->GetPoint (i) << " (grid : " << (mesh->GetPoint (i)->GetGlobalIndex () < mesh->GetNumberOfCartesianPoints ())<< ")" << std::endl;
    //    }

    Vector value (mesh->GetNumberOfTotalPoints ());
    value.setZero ();

    for (int i = 0; i < phi_value.size (); ++i)
    {
        value.coeffRef (i) = phi_value.coeffRef (i);
    }
    //    for (int k = 0; k < mesh->Get_Nz (); ++k)
    //        for (int j = 0; j < mesh->Get_Ny (); ++j)
    //            for (int i = 0; i < mesh->Get_Nx (); ++i)
    //            {

    //                value.coeffRef (i) = mesh->GetPoint (i)->GetGlobalIndex ();
    //            }


//    value = FunToVec (mesh, phi, phi_prim);
    Writer * writer = new Writer (mesh);

    writer->SetFilename ("test");
    writer->SetVectorNumerical (&value);
    writer->SetWriteBothDomainsOn ();
    writer->WriteNow ();

    delete mesh;
    delete writer;

    //    Vector u (mesh->GetNumberOfTotalPoints ());
    //    u.setOnes ();

    //    Writer writer (mesh);
    //    writer.SetVectorAnalytical (&u);

    //    writer.SetFilename ("Test");
    //    writer.WriteNow ();

    return 0;
}

//Point phi (double theta, double t)
//{
//    (void)t;

//    Point p;
//    p.x = 0.02 * std::sqrt (5) + (0.5 + 0.2 * std::sin (5 * theta)) * std::cos (theta);
//    p.y = 0.02 * std::sqrt (5) + (0.5 + 0.2 * std::sin (5 * theta)) * std::sin (theta);
//    p.z = 0;

//    return p;
//}

//Point phi_prim (double theta, double t)
//{
//    (void)t;

//    Point p;
//    p.x = std::cos (5 * theta) * std::cos (theta) - (0.5 + 0.2 * std::sin (5 * theta)) * std::sin (theta);
//    p.y = std::cos (5 * theta) * std::sin (theta) + (0.5 + 0.2 * std::sin (5 * theta)) * std::cos (theta);
//    p.z = 0;

//    return p;
//}

Point phi (double theta, double t)
{
    (void)t;

    Point p;
    p.x = 0.6 * std::cos (theta) - 0.3 * std::cos (3 * theta);
    p.y = 1.5 + 0.7 * std::sin (theta) -  0.07 * std::sin (3 * theta) + 0.2 * std::sin (7 * theta);
    p.z = 0;

    return p;
}

Point phi_prim (double theta, double t)
{
    (void)t;

    Point p;
    p.x =  - 0.6 * std::sin (theta) + 0.9 * std::sin(3 * theta);
    p.y = 0.7 * std::cos (theta) -  0.21 * std::cos (3 * theta) + 1.4 * std::cos (7 * theta);
    p.z = 0;

    return p;
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
