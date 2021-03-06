#include "../O2FID/O2FID.h"

#include <iostream>
#include <cmath>
#include "Eigen/SparseCore"

#define SPACER std::left << std::setw(20)


double phi (Point p, double t = 0); // fonction levelset
double f (Point a, double t = 0.); // fonction de second membre
double u (Point a, double t = 0.);

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "            EXAMPLE 1 - O2FID            " << std::endl;
    std::cout << "-----------------------------------------" << std::endl;

    std::vector<int> listNx = {41, 81, 161}; // liste des Nx
//    std::vector<int> listNx = {10, 20, 40, 80, 160, 320, 640, 1280}; // liste des Nx

    std::vector<double> err_l1 = {}; // erreur l1
    std::vector<double> err_linf = {}; // erreur linf
    std::vector<double> err_rela = {}; // erreur relative
    std::vector<double> h = {}; // pas d'espace h

    for (size_t idx = 0; idx < listNx.size (); ++idx)
    {
        int Nx = listNx.at (idx);

        // Construction du MESH
        Mesh* mesh = new Mesh ();

        mesh->SetBounds (new Point(-0.5, 0, 0), new Point(0.5, 0, 0));
        mesh->Set_Nx(Nx);
//        mesh->Set_Ny(Ny);
//        mesh->Set_Nz(Nz);
        mesh->Build ();

        // Construction de vecteur phi fonction de levelset
        Vector phi_vec = FunToVec (mesh, phi);

        // Ajout des points de bord
        std::vector<int> listPoint = MakeBorderPoints (mesh, &phi_vec);

        mesh->Print ();

        // Construction de la matrice du Laplacien
        Matrix A = Laplacian (mesh);

        // Construction du vecteur de second membre
        Vector b = FunToVec (mesh, f);

        // Imposition de Dirichlet sur les points situés sur le bord du domaine Omega
        ImposeDirichlet (mesh, &A, &b, u, listPoint);

        // Imposition de Dirichlet sur les points situés sur le bord du grand domaine
        //
        //  ImposeDirichlet (mesh, &A, &b, u, {0, mesh->GetNumberOfCartesianPoints ()-1});
        //

        // Vecteur de solution numérique
        Vector u_num = Solve (A, b, IMPLICIT);

        // Vecteur de solution analytique
        Vector u_ana = FunToVec (mesh, u);

        // Transforme les deux vecteurs en 0 sur EXTERNOMEGA pour calculer des erreurs
        mesh->MakeZeroOnExternOmegaInVector (&u_ana);
        mesh->MakeZeroOnExternOmegaInVector (&u_num);

        // Erreur en valeur absolue
        Vector err_abs = GetErrorAbs (mesh, u_ana, u_num);

        // Erreur l1
        err_l1.push_back (GetErrorl1 (mesh, u_ana, u_num));

        // Erreur linf
        err_linf.push_back (GetErrorlinf (mesh, u_ana, u_num));

        // Erreur relative
        err_rela.push_back (GetErrorRela (mesh, u_ana, u_num));

        // Pas h du maillage (rayon de la boule ?)
        Point p = {mesh->Get_hx (), mesh->Get_hy (), mesh->Get_hz ()};
        h.push_back (std::sqrt(p|p));

        // Écriture dans des fichiers
        Writer writer (mesh);
        writer.SetFilename (std::string ("example_1_") + std::to_string (Nx));
//        writer.SetCurrentIteration (0); // Itérations lorsqu'il y a du temps
        writer.SetVectorNumerical (&u_num);
        writer.SetVectorAnalytical (&u_ana);
        writer.SetWriteBothDomainsOn (); // Écrire sur le domaine entier ?
        writer.SetVectorErrorAbs (&err_abs);

        writer.WriteNow ();

        delete mesh;
    }

    std::cout << std::endl;

    std::cout << "#Summary " << std::endl;

    std::cout << "Nx            : " << listNx << std::endl;
    std::cout << "h            : " << h << std::endl;
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

    return fabs(p.x) - 0.313;
}

double f (Point a, double t)
{
    (void)t;

    double c = std::cos(2. * M_PI * a.x);
    double s = std::sin(2. * M_PI * a.x);

    return 8. * (1. - 2. * M_PI * M_PI * a.x * a.x) * s + 32. * M_PI * a. x * c;
}

double u (Point a, double t)
{
    (void)t;
    return 4. * a.x * a.x * std::sin(2. * M_PI * a.x);
}
