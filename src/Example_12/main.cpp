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


    int Nx = 200;
    int Ny = Nx;
    //    int Nz = Nx;

    // Construction du MESH
    Mesh* mesh = new Mesh ();

    mesh->SetBounds (new Point(-1., -1., 0), new Point(1., 1., 0));
    //    mesh->SetBounds (new Point(-0.2, -0.2, -0.2), new Point(0.2, 0.2, 0.2));

    mesh->Set_Nx(Nx);
    mesh->Set_Ny(Ny);
    //    mesh->Set_Nz(Nz);
    mesh->Build ();

    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hx ();
    double hz = mesh->Get_hz ();

    //        double dt = 0.4 * std::min (1. / (1./ hx + 1./ hy), hx);
    //            double dt = 0.5 / (1./ hx + 1./ hy);
    //    double dt = 0.5 * hx * hx;
    double dt = 0.01 * hx * hx;



    // Construction des vecteurs phi fonction de levelset
    Vector phi_vec = FunToVec (mesh, phi, 0.);

    // Ajout des points de bord
    std::vector<int> listPoint = MakeBorderPoints (mesh, &phi_vec);
    phi_vec = FunToVec (mesh, phi, 0.);

    mesh->Print ();

    // T_num
    Vector T_num = FunToVec (mesh, T, 0.);

    // Field;
    Field* field = GetWField (mesh, &phi_vec, &T_num, &listPoint, 1.);
    Vector NormGradPhi = ComputeNorm(&field->GradPhi);

    // Écriture dans des fichiers
    Writer writer (mesh);
    writer.SetFilename (std::string ("example_12_") + std::to_string (Nx));

    double theta = 1.;
    int IterMax = 50;

    bool mybool = true;

    double time = 1.;

    for (int timer = 0; timer <= IterMax && mybool; ++timer)
    {
        time += dt;

        std::cout << "-----------------------------------------------" << std::endl;
        std::cout << INDENT << "Time : " << time << " [dt = " << dt << "](" << timer << "/" << IterMax << ") --- theta : " << theta << std::endl;
        std::cout << "-----------------------------------------------" << std::endl;

        mesh->Print ();
        if (timer > 0)
        {
            if (int(listPoint.size ()) <= mesh->GetNumberOfCartesianPoints ())
            {
                // Construction de la matrice du Laplacien
                Matrix As = -1. * Laplacian (mesh);
                Matrix At = 0. * As;

                At.pruned ();
                At.setIdentity ();
                At *= 1./ dt;

                Extrapole (mesh, &T_num);

                // Construction du vecteur de second membre
                Vector b = (At + (1. - theta) * As) * T_num;

                // Matrix
                Matrix A = At + theta * As;

                ImposeDirichlet (mesh, &A, &b, T, mesh->GetListOfIndexPoints ());

                T_num = Solve (A, b, IMPLICIT);

            } else
            {
                mybool = false;
            }
        }

        // Extrapole partie
        Extrapole (mesh, &NormGradPhi);
        Extrapole (mesh, &field->W);
        Extrapole (mesh, &field->Normals);
        Extrapole (mesh, &field->GradPhi);
        Extrapole (mesh, &field->GradTemperature);

        // Écriture
        writer.SetCurrentIteration (timer);
        writer.SetVectorNumerical (&T_num);
        writer.SetVectorPhi (&phi_vec);
        writer.SetNormPhi (&NormGradPhi);
        writer.SetVectorW_new (&field->W);
        writer.SetVectorNormals (&field->Normals);
        writer.SetVectorGradPhi (&field->GradPhi);
        writer.SetVectorGradTemperature (&field->GradTemperature);
        writer.WriteNow ();


        // Compute for Next iteration
        delete field;
        listPoint = mesh->GetListOfIndexPoints ();

        field = GetWField (mesh, &phi_vec, &T_num, &listPoint, 1.);

//        dt = Compute_dt (mesh, field);

        NormGradPhi = ComputeNorm(&field->GradPhi);

        phi_vec = IteratePhi (mesh, &field->W, dt, &phi_vec);

        mesh->RemoveAllNotCartesianPoints ();

        listPoint = MakeBorderPoints (mesh, &phi_vec);

        phi_vec = ReInitPhi (mesh, &phi_vec, 3, TVD_RK4);

        Extrapole (mesh, &phi_vec);

    }

    delete field;
    delete mesh;


    return 0;
}

double phi (Point p, double t)
{
    (void)t;

    return std::sqrt(p|p) - 0.05;
    //    return std::sqrt(p|p) - (0.4 + 0.1 * std::sin (5. * std::atan(p.y / fabs(p.x))));

}

double T (Point p, double t)
{
    (void)t;

    switch (p.GetLocate ())
    {
    case ON_DOMAIN_INTERN_OMEGA:
        //        std::cout << "Request Intern" << std::endl;
        return 0.;
    case ON_BORDER_OMEGA:
        //        std::cout << "Request Border" << std::endl;
        return -0.5;
    case ON_DOMAIN_EXTERN_OMEGA:
        //        std::cout << "Request Extern" << std::endl;
        return -0.5;
    }

    //    return -0.5;
}


