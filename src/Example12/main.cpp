#include "../O2FID/O2FID.h"

#include <iostream>
#include <cmath>
#include "Eigen/SparseCore"
#include "Eigen/Dense"
#include "Eigen/Core"


#define SPACER std::left << std::setw(20)


double phi (Point p, double t = 0.); // fonction levelset
double T (Point p, double t = 0.);

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "            EXAMPLE 12 - O2FID            " << std::endl;
    std::cout << "-----------------------------------------" << std::endl;


    int Nx = 400;
    int Ny = Nx;


    // Construction du MESH
    Mesh* mesh = new Mesh ();

    mesh->SetBounds (new Point(-2, -2, 0), new Point(2, 2, 0));
    mesh->Set_Nx(Nx);
    mesh->Set_Ny(Ny);
    //        mesh->Set_Nz(Nz);
    mesh->Build ();

    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hx ();
    double hz = mesh->Get_hz ();
    int N = mesh->GetNumberOfTotalPoints ();
    double dt = 0.5 * std::min (1. / (1./ hx + 1./ hy), hx);
    //    double dt = 0.5 / (1./ hx + 1./ hy);


    // Construction des vecteurs phi fonction de levelset
    Vector phi_vec = FunToVec (mesh, phi);

    // Ajout des points de bord
    std::vector<int> listPoint = MakeBorderPoints (mesh, &phi_vec);
    phi_vec = FunToVec (mesh, phi);

    mesh->Print ();

    // Construction de la matrice du Laplacien
    Matrix As;
    Matrix At;

    // Construction du vecteur de second membre
    Vector b;

    // Matrix
    Matrix A;

    // T_num
    Vector T_num = FunToVec (mesh, T);

    // W
    std::vector<Point*> Wold, Wnew;
    Wnew = GetWField (mesh, &phi_vec, &T_num, &listPoint, 1.);
    Wold = Wnew;

    // Écriture dans des fichiers
    Writer writer (mesh);
    writer.SetFilename (std::string ("example_12_") + std::to_string (Nx));
    writer.SetCurrentIteration (0); // Itérations lorsqu'il y a du temps
    writer.SetVectorNumerical (&T_num);
    writer.SetVectorW_new (&Wnew);
    writer.SetWriteBothDomainsOn (); // Écrire sur le domaine entier ?
    writer.SetVectorPhi (&phi_vec);

    writer.WriteNow ();



    double theta = 1.;
    int IterMax = 50;

    bool mybool = true;

    for (int timer = 1; timer <= IterMax && mybool; ++timer)
    {
        double time = timer * dt;

        std::cout << "-----------------------------------------------" << std::endl;
        std::cout << INDENT << "Time : " << time << " (" << timer << "/" << IterMax << ") --- theta : " << theta << std::endl;
        std::cout << "-----------------------------------------------" << std::endl;

        phi_vec = IteratePhi (mesh, &Wold, dt, &phi_vec);
//        ReInitPhi (mesh, &phi_vec, &listPoint, dt);

        mesh->RemoveAllNotCartesianPoints ();

        listPoint = MakeBorderPoints (mesh, &phi_vec);

        mesh->Print ();

        N = mesh->GetNumberOfTotalPoints ();

        Vector Temp = T_num;

        if (int(listPoint.size ()) <= mesh->GetNumberOfCartesianPoints ())
        {
            As = Laplacian (mesh);
            At.resize (N, N);
            At.setIdentity ();
            At *= 1./ dt;

            Extrapole (mesh, &Temp);

            //        A_cn = (0.5 * As + At);
            //        b = At * T_num;
            //        b -= 0.5 * T_num;

            b = (At + (1. - theta) * As) * T_num;
            A = At - theta * As;

            ImposeDirichlet (mesh, &A, &b, T, listPoint);


            T_num = Solve (A, b, IMPLICIT);
        } else
        {
            mybool = false;
        }


        Wold = Wnew;
        Extrapole (mesh, &Wold);

        Wnew = GetWField (mesh, &phi_vec, &T_num, &listPoint, 1.);

        writer.SetCurrentIteration (timer);
        writer.SetVectorNumerical (&T_num);
        writer.SetVectorAnalytical (&Temp);
        writer.SetVectorPhi (&phi_vec);
        writer.SetVectorW_new (&Wnew);
        writer.SetVectorW_old (&Wold);

        writer.WriteNow ();
    }

    delete mesh;


    return 0;
}

double phi (Point p, double t)
{
    (void)t;

    return (p|p) - 0.4*0.4;
}

double T (Point p, double t)
{
    (void)t;

    switch (p.GetLocate ())
    {
    case ON_DOMAIN_INTERN_OMEGA:
        //        std::cout << "Request Intern" << std::endl;
        return 0;
    case ON_BORDER_OMEGA:
        //        std::cout << "Request Border" << std::endl;
        return -0.5;
    case ON_DOMAIN_EXTERN_OMEGA:
        //        std::cout << "Request Extern" << std::endl;
        return -0.5;
    }

    return -0.5;
}


