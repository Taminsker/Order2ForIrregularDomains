#include "../O2FID/O2FID.h"

#include <iostream>
#include <cmath>
#include "Eigen/SparseCore"

#define SPACER std::left << std::setw(20)


double phi (Point p, double t = 0); // fonction levelset
double T (Point a, double t = 0.);

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "           EXAMPLE 10 - O2FID            " << std::endl;
    std::cout << "-----------------------------------------" << std::endl;

    std::vector<int> listNx = {41, 81, 161}; // liste des Nx
    std::vector<double> err_l1 = {}; // erreur l1
    std::vector<double> err_linf = {}; // erreur linf
    std::vector<double> err_rela = {}; // erreur relative
    std::vector<double> h = {}; // pas d'espace h

    int Nt = 30;

    for (size_t idx = 0; idx < listNx.size (); ++idx)
    {
        int Nx = listNx.at (idx);

        // Construction du MESH
        Mesh* mesh = new Mesh ();

        mesh->SetBounds (new Point(), new Point(1));
        mesh->Set_Nx(Nx);

        mesh->Build ();

        double dt = mesh->Get_hx ();

        // Construction de vecteur phi fonction de levelset
        Vector phi_vec = FunToVec (mesh, phi);

        // Ajout des points de bord
        std::vector<int> listPoint = MakeBorderPoints (mesh, &phi_vec);
        //        ExtendToNeighbours (mesh, &listPoint);

        mesh->Print ();

        // Vecteur de solution
        Vector u_num, u_ana;
        u_num = u_ana = FunToVec (mesh, T, 0.);

        // Matrice
        int N = mesh->GetNumberOfTotalPoints ();

        Matrix A = -1. * Laplacian (mesh);
        //        RemovePeriodicity (mesh, &A);

        Matrix Id (N, N);
        Id.setIdentity ();

        Matrix Al = (1./ dt) * Id + A;
        Matrix Ar = (1./ dt) * Id;

        Vector b;

        // Field;
        Field field;
        std::vector<Point*> Wold;
        GetWField (&field, mesh, &phi_vec, &u_num, &listPoint, 1.);
        Wold = field.W;

        Vector err_abs = GetErrorAbs (mesh, u_ana, u_num);

        // Écriture dans des fichiers
        Writer writer (mesh);
        writer.SetFilename (std::string ("example_10_") + std::to_string (Nx));
        writer.SetCurrentIteration (0); // Itérations lorsqu'il y a du temps
        writer.SetVectorNumerical (&u_num);
        writer.SetVectorW_new (&field.W);

        writer.SetWriteBothDomainsOn (); // Écrire sur le domaine entier ?
        writer.SetVectorPhi (&phi_vec);
        writer.SetVectorNormals (&field.Normals);
        writer.SetVectorGradPhi (&field.GradPhi);
        writer.SetVectorGradTemperature (&field.GradTemperature);
        writer.SetVectorErrorAbs (&err_abs);

        writer.WriteNow ();

        for (int it = 1; it <= Nt; ++it)
        {
            double t = double(it) * dt;

            std::cout << "-----------------------------------------------" << std::endl;
            std::cout << INDENT << "Time : " << t << " (" << it << "/" << Nt << ") " << std::endl;
            std::cout << "-----------------------------------------------" << std::endl;

            phi_vec = IteratePhi (mesh, &Wold, dt, &phi_vec);
            //        ReInitPhi (mesh, &phi_vec, &listPoint, dt);

            mesh->RemoveAllNotCartesianPoints ();

            listPoint = MakeBorderPoints (mesh, &phi_vec);

            Matrix As = Al;
            Vector b = Ar * u_num;

            ImposeDirichlet (mesh, &As, &b, T, listPoint, t);

            u_num = Solve (As, b, IMPLICIT);

            Extrapole (mesh, &u_num);
            for (int idx : listPoint)
                u_num.coeffRef (idx) = T(*mesh->GetPoint (idx), t);

            u_ana = FunToVec (mesh, T, t);

            // Erreur en valeur absolue
            err_abs = GetErrorAbs (mesh, u_ana, u_num);

            if (it == Nt - 1)
            {
                // Erreur l1
                err_l1.push_back (GetErrorl1 (mesh, u_ana, u_num));

                // Erreur linf
                err_linf.push_back (GetErrorlinf (mesh, u_ana, u_num));

                // Erreur relative
                err_rela.push_back (GetErrorRela (mesh, u_ana, u_num));

                // Pas h du maillage (rayon de la boule ?)
                Point p = {mesh->Get_hx (), mesh->Get_hy (), mesh->Get_hz ()};
                h.push_back (std::sqrt(p|p));
            }

            // Écriture dans des fichiers

            Wold = field.W;
            Extrapole (mesh, &Wold);

    //        Wnew = GetWField (mesh, &phi_vec, &T_num, &listPoint, 1.);


            Extrapole (mesh, &phi_vec);
            writer.SetCurrentIteration (it);
            writer.SetVectorNumerical (&u_num);
    //        writer.SetVectorAnalytical (&Temp);
            writer.SetVectorPhi (&phi_vec);
//            writer.SetNormPhi (&GradPhi);
            writer.SetVectorW_new (&field.W);
            writer.SetVectorNormals (&field.Normals);
            writer.SetVectorGradPhi (&field.GradPhi);
            writer.SetVectorGradTemperature (&field.GradTemperature);
            writer.SetVectorW_old (&Wold);
            writer.WriteNow ();
        }

        delete mesh;
    }

    std::cout << std::endl;

    std::cout << "#Summary " << std::endl;

    std::cout << "Nx            : " << listNx << std::endl;
    std::cout << "l1-error      : " << err_l1 << std::endl;
    std::cout << "Order         : " << Order(err_l1, h) << std::endl;
    std::cout << "linf-error    : " << err_linf << std::endl;
    std::cout << "Order         : " << Order(err_linf, h) << std::endl;
    std::cout << "rela-error    : " << err_rela << std::endl;
    std::cout << "Order         : " << Order(err_rela, h) << std::endl;

    return 0;
}

double phi (Point p, double t)
{
    (void)t;

    return p.x - 0.5;
}

double T (Point a, double t)
{
    (void)t;

    if (a.GetLocate () == ON_BORDER_OMEGA || a.GetLocate () == ON_DOMAIN_INTERN_OMEGA)
        return std::exp(t - a.x + 0.5) - 1.;
    else
        return 0.;
}
