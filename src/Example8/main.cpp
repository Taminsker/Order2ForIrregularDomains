#include "../O2FID/O2FID.h"

#include <iostream>
#include <cmath>
#include "Eigen/SparseCore"

#define SPACER std::left << std::setw(20)


Point phi (double theta, double t = 0); // fonction levelset
Point phigrad (double theta, double t=0);
double f (Point a, double t = 0.); // fonction de second membre
double u (Point a, double t = 0.);

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "            EXAMPLE 8 - O2FID            " << std::endl;
    std::cout << "-----------------------------------------" << std::endl;

    std::vector<int> listNx = {81, 161, 321}; // liste des Nx
    std::vector<double> err_l1 = {}; // erreur l1
    std::vector<double> err_linf = {}; // erreur linf
    std::vector<double> err_rela = {}; // erreur relative
    std::vector<double> h = {}; // pas d'espace h

    for (size_t idx = 0; idx < listNx.size (); ++idx)
    {
        int Nx = listNx.at (idx);
        int Ny = Nx;
        int Nt = 20;

        // Construction du MESH
        Mesh* mesh = new Mesh ();

        //pas de temps
        double dt = mesh->Get_hx();
        double dt_BE2 = mesh->Get_hx() * mesh->Get_hx();
        int N = mesh->GetNumberOfTotalPoints();

        mesh->SetBounds (new Point(-1.0, -1.0, 0), new Point(1.0, 1.0, 0));
        mesh->Set_Nx(Nx);
        mesh->Set_Ny(Ny);
        mesh->Build ();

        // Construction de vecteur phi fonction de levelset
        Vector phi_vec = FunToVec (mesh, phi, phigrad);

        // Ajout des points de bord
        std::vector<int> listPoint = MakeBorderPoints (mesh, &phi_vec);

        mesh->Print ();

        // Construction de la matrice du Laplacien (Crank Nicolson)
        Matrix Id (N, N); // création de la matrice identité de taille N x N
        Id.setIdentity();
        Matrix A_CK = (1/dt) * Id + (1/2) * Laplacian (mesh);

        // Construction de la matrice du Laplacien (Backward-Euler)
        Matrix A_BE = Laplacian (mesh);

        // Construction du vecteur de second membre
        Vector b = FunToVec (mesh, f);

        // Imposition de Dirichlet sur les points situés sur le bord du domaine Omega
        ImposeDirichlet (mesh, &A_CK, &b, u, listPoint);
        ImposeDirichlet (mesh, &A_BE, &b, u, listPoint);

        // Imposition de Dirichlet sur les points situés sur le bord du grand domaine
        //
        //  ImposeDirichlet (mesh, &A, &b, u, {0, mesh->GetNumberOfCartesianPoints ()-1});
        //

        // CRANK-NICOLSON
        for (int idt = 0; idt < Nt; idt++) // boucle sur les points temporelles
        {
            // Vecteur de solution numérique
            Vector u_num_CK = Solve (A_CK, b, IMPLICIT);
            // U_num deviens b pour le temps suivant
            Vector b = u_num_CK;

            // Vecteur de solution analytique à t'instant t = i*dt
            Vector u_ana_CK = FunToVec (mesh, u, idt*dt);

            // Transforme les deux vecteurs en 0 sur EXTERNOMEGA pour calculer des erreurs
            mesh->MakeZeroOnExternOmegaInVector (&u_ana_CK);
            mesh->MakeZeroOnExternOmegaInVector (&u_num_CK);

            // Erreur en valeur absolue
            Vector err_abs_CK = GetErrorAbs (mesh, u_ana_CK, u_num_CK);

            // Erreur l1
            err_l1.push_back (GetErrorl1 (mesh, u_ana_CK, u_num_CK));

            // Erreur linf
            err_linf.push_back (GetErrorlinf (mesh, u_ana_CK, u_num_CK));

            // Erreur relative
            err_rela.push_back (GetErrorRela (mesh, u_ana_CK, u_num_CK));

            // Pas h du maillage (rayon de la boule ?)
            Point p = {mesh->Get_hx (), mesh->Get_hy (), mesh->Get_hz ()};
            h.push_back (std::sqrt(p|p));

            // Écriture dans des fichiers
            Writer writer (mesh);
            writer.SetFilename (std::string ("example_8_Crank_Nicolson_Nx_") + std::to_string (Nx) + std::string ("_dt_") + std::to_string (idt*dt));
            writer.SetCurrentIteration (idt); // Itérations lorsqu'il y a du temps
            writer.SetVectorNumerical (&u_num_CK);
            writer.SetVectorAnalytical (&u_ana_CK);
            writer.SetWriteBothDomainsOn (); // Écrire sur le domaine entier ?
            writer.SetVectorErrorAbs (&err_abs_CK);

            writer.WriteNow ();
        }

        std::cout << std::endl;

        std::cout << "#Summary Crank-Nicolson " << std::endl;

        std::cout << "Nx            : " << listNx << std::endl;
        std::cout << "l1-error      : " << err_l1 << std::endl;
        std::cout << "Order         : " << Order(err_l1, h) << std::endl;
        std::cout << "linf-error    : " << err_linf << std::endl;
        std::cout << "Order         : " << Order(err_linf, h) << std::endl;
        std::cout << "rela-error    : " << err_rela << std::endl;
        std::cout << "Order         : " << Order(err_rela, h) << std::endl;

        // BACKWARD-EULER 1 dt ~ dx
        for (int idt = 0; idt < Nt; idt++) // boucle sur les points temporelles
        {
            // Vecteur de solution numérique
            Vector u_num_BE1 = Solve (A_BE, b, IMPLICIT);
            // U_num deviens b pour le temps suivant
            Vector b = u_num_BE1;

            // Vecteur de solution analytique à t'instant t = i*dt
            Vector u_ana_BE1 = FunToVec (mesh, u, idt*dt);

            // Transforme les deux vecteurs en 0 sur EXTERNOMEGA pour calculer des erreurs
            mesh->MakeZeroOnExternOmegaInVector (&u_ana_BE1);
            mesh->MakeZeroOnExternOmegaInVector (&u_num_BE1);

            // Erreur en valeur absolue
            Vector err_abs_BE1 = GetErrorAbs (mesh, u_ana_BE1, u_num_BE1);

            // Erreur l1
            err_l1.push_back (GetErrorl1 (mesh, u_ana_BE1, u_num_BE1));

            // Erreur linf
            err_linf.push_back (GetErrorlinf (mesh, u_ana_BE1, u_num_BE1));

            // Erreur relative
            err_rela.push_back (GetErrorRela (mesh, u_ana_BE1, u_num_BE1));

            // Pas h du maillage (rayon de la boule ?)
            Point p = {mesh->Get_hx (), mesh->Get_hy (), mesh->Get_hz ()};
            h.push_back (std::sqrt(p|p));

            // Écriture dans des fichiers
            Writer writer (mesh);
            writer.SetFilename (std::string ("example_8_Backward_Euler1_Nx_") + std::to_string (Nx) + std::string ("_dt_") + std::to_string (idt*dt));
            writer.SetCurrentIteration (idt); // Itérations lorsqu'il y a du temps
            writer.SetVectorNumerical (&u_num_BE1);
            writer.SetVectorAnalytical (&u_ana_BE1);
            writer.SetWriteBothDomainsOn (); // Écrire sur le domaine entier ?
            writer.SetVectorErrorAbs (&err_abs_BE1);

            writer.WriteNow ();
        }

        std::cout << std::endl;

        std::cout << "#Summary Backward-Euler dt ~ dx " << std::endl;

        std::cout << "Nx            : " << listNx << std::endl;
        std::cout << "l1-error      : " << err_l1 << std::endl;
        std::cout << "Order         : " << Order(err_l1, h) << std::endl;
        std::cout << "linf-error    : " << err_linf << std::endl;
        std::cout << "Order         : " << Order(err_linf, h) << std::endl;
        std::cout << "rela-error    : " << err_rela << std::endl;
        std::cout << "Order         : " << Order(err_rela, h) << std::endl;

        // BACKWARD-EULER 2 dt ~ dx*dx
        for (int idt = 0; idt < Nt; idt++) // boucle sur les points temporelles
        {
            // Vecteur de solution numérique
            Vector u_num_BE2 = Solve (A_BE, b, IMPLICIT);
            // U_num deviens b pour le temps suivant
            Vector b = u_num_BE2;

            // Vecteur de solution analytique à t'instant t = i*dt
            Vector u_ana_BE2 = FunToVec (mesh, u, idt*dt_BE2);

            // Transforme les deux vecteurs en 0 sur EXTERNOMEGA pour calculer des erreurs
            mesh->MakeZeroOnExternOmegaInVector (&u_ana_BE2);
            mesh->MakeZeroOnExternOmegaInVector (&u_num_BE2);

            // Erreur en valeur absolue
            Vector err_abs_BE2 = GetErrorAbs (mesh, u_ana_BE2, u_num_BE2);

            // Erreur l1
            err_l1.push_back (GetErrorl1 (mesh, u_ana_BE2, u_num_BE2));

            // Erreur linf
            err_linf.push_back (GetErrorlinf (mesh, u_ana_BE2, u_num_BE2));

            // Erreur relative
            err_rela.push_back (GetErrorRela (mesh, u_ana_BE2, u_num_BE2));

            // Pas h du maillage (rayon de la boule ?)
            Point p = {mesh->Get_hx (), mesh->Get_hy (), mesh->Get_hz ()};
            h.push_back (std::sqrt(p|p));

            // Écriture dans des fichiers
            Writer writer (mesh);
            writer.SetFilename (std::string ("example_8_Backward_Euler2_Nx_") + std::to_string (Nx) + std::string ("_dt_") + std::to_string (idt*dt_BE2));
            writer.SetCurrentIteration (idt); // Itérations lorsqu'il y a du temps
            writer.SetVectorNumerical (&u_num_BE2);
            writer.SetVectorAnalytical (&u_ana_BE2);
            writer.SetWriteBothDomainsOn (); // Écrire sur le domaine entier ?
            writer.SetVectorErrorAbs (&err_abs_BE2);

            writer.WriteNow ();
        }

        delete mesh;
    }

    std::cout << std::endl;

    std::cout << "#Summary Backward-Euler dt ~ dx^2 " << std::endl;

    std::cout << "Nx            : " << listNx << std::endl;
    std::cout << "l1-error      : " << err_l1 << std::endl;
    std::cout << "Order         : " << Order(err_l1, h) << std::endl;
    std::cout << "linf-error    : " << err_linf << std::endl;
    std::cout << "Order         : " << Order(err_linf, h) << std::endl;
    std::cout << "rela-error    : " << err_rela << std::endl;
    std::cout << "Order         : " << Order(err_rela, h) << std::endl;

    return 0;
}

Point phi (double theta, double t)
{
    (void)t;

    Point p;
    p.x = 0.02 * std::sqrt (5.) + (0.5 + 0.2 * std::sin(5. * theta)) * std::cos (theta);
    p.y = 0.02 * std::sqrt (5.) + (0.5 + 0.2 * std::sin(5. * theta)) * std::sin (theta);
    p.z = 0.;

    return p;
}

Point phigrad (double theta, double t)
{
    (void)t;

    Point p;
    p.x = 0.2 * std::cos(5*theta) * std::cos(theta) - (0.5 + 0.2 * std::sin(5*theta)) * std::sin(theta);
    p.y = 0.2 * std::cos(5*theta) * std::sin(theta) + (0.5 + 0.2 * std::sin(5*theta)) * std::cos(theta);
    p.z = 0.;

    return p;
}

double f (Point a, double t)
{
    return 0;
}

double u (Point a, double t)
{
    return std::exp(-2 * t) * std::sin(a.x) * std::sin(a.y);
}
