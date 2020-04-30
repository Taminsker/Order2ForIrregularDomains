#include "../O2FID/O2FID.h"

#include <iostream>
#include <cmath>
#include "Eigen/SparseCore"

#define SPACER std::left << std::setw(20)


double phi (Point p, double t = 0); // fonction levelset
double beta (Point a, double t = 0);
double f (Point a, double t = 0.); // fonction de second membre
double u (Point a, double t = 0.);

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "         EXAMPLE TVD WENO - O2FID        " << std::endl;
    std::cout << "-----------------------------------------" << std::endl;

//    std::vector<int> listNx = {21, 41, 81, 201}; // liste des Nx
    std::vector<int> listNx = {41}; // liste des Nx

    std::vector<double> err_l1 = {}; // erreur l1
    std::vector<double> err_linf = {}; // erreur linf
    std::vector<double> err_rela = {}; // erreur relative
    std::vector<double> h = {}; // pas d'espace h

    for (size_t idx = 0; idx < listNx.size (); ++idx)
    {
        int Nx = listNx.at (idx);
        int Ny = Nx;
        int Nz = 1;
        int Nt = 20;

        // Construction du MESH
        Mesh* mesh = new Mesh ();

        mesh->SetBounds (new Point(-1., -1.), new Point(1, 1));
        mesh->Set_Nx(Nx);
        mesh->Set_Ny(Ny);
        mesh->Build ();

        double dt = 0.3 * mesh->Get_hx ();

        // Construction de vecteur phi fonction de levelset
        Vector phi_vec = FunToVec (mesh, phi);

        // Ajout des points de bord
        std::vector<int> listPoint = MakeBorderPoints (mesh, &phi_vec);
        //        ExtendToNeighbours (mesh, &listPoint);

        mesh->Print ();

        //        std::cout << Gradient (mesh, ORDER_2_CENTRAL) << std::endl;

        //        delete mesh;
        //        return 0;

        // Construction du vecteur de second membre
        std::vector<Point*> W;
        int N = mesh->GetNumberOfTotalPoints ();
        W.resize (size_t(N));

        for (size_t i = 0; i < size_t(N); ++i)
            W.at (i) = new Point(M_PI/2., -0.90);

        // Vecteur de solution
        Vector u_num, u_ana;
        u_num = u_ana = FunToVec (mesh, u, 0.);

        Vector err_abs = GetErrorAbs (mesh, u_ana, u_num);

        // Écriture dans des fichiers
        Writer writer (mesh);
        writer.SetFilename (std::string ("example_TVD_WENO_") + std::to_string (Nx));
        writer.SetCurrentIteration (0); // Itérations lorsqu'il y a du temps
        writer.SetVectorNumerical (&u_num);
        writer.SetVectorAnalytical (&u_ana);
        writer.SetWriteBothDomainsOn (); // Écrire sur le domaine entier ?
        writer.SetVectorErrorAbs (&err_abs);

        writer.WriteNow ();

        Matrix A = Gradient (mesh, ORDER_2_FORWARD);

        for (int it = 1; it <= Nt; ++it)
        {
            double t = double(it) * dt;

            std::cout << "t : " << t << std::endl;

            if (true)
            {
                // TVD scheme
                Vector u_n = u_num;
                Vector u_1 = u_n + dt * Weno (mesh, &u_n, &W);
                Vector u_2 = 3./4. * u_n + 1./ 4. * u_1 + 1./ 4. * dt * Weno (mesh, &u_1, &W);
                u_num = 1./ 3. * u_n + 2. / 3. * u_2 + 2./ 3. * dt * Weno (mesh, &u_2, &W);
            }

            if (false)
            {
                u_num = u_num + dt * Weno(mesh, &u_num, &W);
                //                u_num = u_num - dt * A;
            }

            if (false)
            {
                Vector tmp = u_num;
                //u_ij^n+1 = u_ij^n + dt * 1/dx(u_i+1^n - u_ij^n)

                for (int i = 0; i < Nx; ++i)
                    for (int j = 0; j < Ny; ++j)
                    {
                        int idx = mesh->GetGlobalIndexOfPoint (i, j);

                        int idx_i = mesh->GetGlobalIndexOfPoint (i+1, j);
                        int idx_j = mesh->GetGlobalIndexOfPoint (i, j+1);

                        u_num.coeffRef (idx) = tmp.coeff (idx) + dt * (1. / mesh->Get_hx () * (tmp.coeff (idx_i) - tmp.coeff (idx)) + 1./mesh->Get_hy () * (tmp.coeff (idx_j) - tmp.coeff (idx)));
                    }

//                u_num = tmp;
            }

            //            Extrapole (mesh, &u_num);
            //            for (int idx : listPoint)
            //                u_num.coeffRef (idx) = u(*mesh->GetPoint (idx), t+dt);

            //            u_ana = FunToVec (mesh, u, t);

            //            mesh->MakeZeroOnExternOmegaInVector (&u_num);
            //            mesh->MakeZeroOnExternOmegaInVector (&u_ana);


            // Erreur en valeur absolue
            //            err_abs = GetErrorAbs (mesh, u_ana, u_num);

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

            writer.SetCurrentIteration (it); // Itérations lorsqu'il y a du temps
            writer.SetVectorNumerical (&u_num);
            //            writer.SetVectorAnalytical (&u_ana);
            writer.SetWriteBothDomainsOn (); // Écrire sur le domaine entier ?
            //            writer.SetVectorErrorAbs (&err_abs);

            writer.WriteNow ();
        }

        AutoClearVector(W);

        delete mesh;
    }

    std::cout << std::endl;

    std::cout << "#Summary " << std::endl;

    std::cout << "Nx            : " << listNx << std::endl;
    std::cout << "Ny            : " << listNx << std::endl;
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
    (void)p;

    //    return EuclidianDist (p, Point()) - 0.5*0.5;

    double eps = 1e-10;

    //    return p.x - .3333;

    if (fabs(p.x - 1.) < eps || fabs(p.y - 1.) < eps || fabs(p.x + 1.) < eps || fabs(p.y + 1.) < eps)
        return 0.;
    else
        return -1.;
}

double f (Point a, double t)
{
    (void)t;
    (void)a;

    return 0.;
}

double u (Point a, double t)
{
    (void)t;

    double eps = 1e-10;

    if (fabs(a.x) < 0.5 && fabs(a.y) < 0.5)
        return std::sin(2. * M_PI * (a.x + t)) * std::sin(2. * M_PI * (a.y + t));
    return 0.;
}
