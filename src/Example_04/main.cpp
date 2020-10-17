#include "../O2FID/O2FID.h"

#include <iostream>
#include <cmath>
#include "Eigen/SparseCore"

#define SPACER std::left << std::setw(20)


Point phi (double theta, double t = 0); // fonction levelset
Point phigrad (double theta, double t = 0);
double beta (Point a, double t);
double f (Point a, double t = 0.); // fonction de second membre
double u (Point a, double t = 0.);

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "            EXAMPLE 4 - O2FID            " << std::endl;
    std::cout << "-----------------------------------------" << std::endl;

    std::vector<int> listNx = {81, 161, 321}; // liste des Nx
    std::vector<int> listNy = {121, 241, 481}; // liste des Nx

//    std::vector<int> listNx = {20, 40, 80, 160, 320, 640}; // liste des Nx
//    std::vector<int> listNy = {30, 60, 120, 240, 480, 960}; // liste des Nx


    std::vector<double> err_l1 = {}; // erreur l1
    std::vector<double> err_linf = {}; // erreur linf
    std::vector<double> err_rela = {}; // erreur relative
    std::vector<double> h = {}; // pas d'espace h

    for (size_t idx = 0; idx < listNx.size (); ++idx)
    {
        int Nx = listNx.at (idx);
        int Ny = listNy.at (idx);

        // Construction du MESH
        Mesh* mesh = new Mesh ();

        mesh->SetBounds (new Point(-1., 0, 0), new Point(1., 3., 0));
        mesh->Set_Nx(Nx);
        mesh->Set_Ny(Ny);
        //        mesh->Set_Nz(Nz);
        mesh->Build ();

        // Construction de vecteur phi fonction de levelset
        Vector phi_vec = FunToVec (mesh, phi, phigrad);

        // Ajout des points de bord
        std::vector<int> listPoint = MakeBorderPoints (mesh, &phi_vec);

        mesh->Print ();

        // Construction de la matrice du Laplacien
        Matrix A = Laplacian (mesh);

        // BETA
        Vector beta_vec = FunToVec (mesh, beta);
        InsertBeta (mesh, &A, &beta_vec);

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
        writer.SetFilename (std::string ("example_4_") + std::to_string (Nx));
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
    std::cout << "Ny            : " << listNy << std::endl;
    std::cout << "h             : " << h << std::endl;
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
    p.x = 0.6 * std::cos (theta) - 0.3 * std::cos (3. * theta);
    p.y = 1.5 + 0.7 * std::sin (theta) - 0.07 * std::sin (3. * theta) + 0.2 * std::sin (7. * theta);
    p.z = 0.;

    return p;
}

Point phigrad (double theta, double t)
{
    (void)t;

    Point p;
    p.x = -0.6 * std::sin (theta) + 0.9 * std::sin (3. * theta);
    p.y = 0.7 * std::cos (theta) - 0.21 * std::cos (3. * theta) + 1.4 * std::cos (7. * theta);
    p.z = 0.;

    return p;
}

double beta (Point a, double t)
{
    (void)t;

    return 2. + std::sin(a.x * a.y);
}

double f (Point a, double t)
{
    (void)t;

    double c_xy = std::cos (a.x * a.y);
    double s_xy = std::sin (a.x * a.y);
    double c_y =  std::cos (a.y);
    double s_y =  std::sin (a.y);

    double y2 = a.y * a.y;
    double x2 = a.x * a.x;

    double x = a.x;
    double y = a.y;

    double value = 0.;

    double e = std::exp(x);

    value = e * (s_xy + 2.) * (2. * s_y + x2 * s_y + 4. * x * s_y + y2);
    value += -1. * e * (s_xy + 2.) * (s_y * x2 - 2.) + y * e * c_xy * (x2 * s_y + 2. * x * s_y + y2);
    value += x * e * c_xy * (c_y * x2 + 2. * y);

    return value;
}

double u (Point a, double t)
{
    (void)t;
    return std::exp (a.x) * (a.x * a.x * std::sin(a.y) + a.y * a.y);
}
