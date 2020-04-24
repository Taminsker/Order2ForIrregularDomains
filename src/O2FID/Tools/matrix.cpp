#include "matrix.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <iostream>


Matrix Laplacian (Mesh * mesh)
{

    std::cout << "# Building the matrix of laplacian" << std::endl;

    // On récupère le nombre de points et les pas d'espace

    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    int Ngrid = mesh->GetNumberOfCartesianPoints (); // Nombre de points de grille uniquement
    int N = mesh->GetNumberOfTotalPoints (); // Points de grille + points rajoutés

    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hy ();
    double hz = mesh->Get_hz ();

    double a = 0., b = 0., c = 0., d = 0.;

    // On définit les coefficients pour implicite par défaut

    if (hx > 1e-10) {b = 1. / (hx * hx);}
    if (hy > 1e-10) {c = 1. / (hy * hy);}
    if (hz > 1e-10) {d = 1. / (hz * hz);}
    a = - 2. * b - 2. * c - 2. * d;

    DIM dim = mesh->GetDimension ();

    Matrix A (N, N);
    A.reserve (Eigen::VectorXi::Constant(N, 15));

    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
            {

                int idx = mesh->GetGlobalIndexOfPoint (i, j, k);

                A.coeffRef (idx, idx) = a;

                int idx_a = mesh->GetGlobalIndexOfPoint (i+1, j, k);
                int idx_b = mesh->GetGlobalIndexOfPoint (i-1, j, k);

                A.coeffRef (idx, idx_a) = b;
                A.coeffRef (idx, idx_b) = b;


                if (dim == DIM_2D || dim == DIM_3D)
                {
                    idx_a = mesh->GetGlobalIndexOfPoint (i, j+1, k);
                    idx_b = mesh->GetGlobalIndexOfPoint (i, j-1, k);

                    A.coeffRef (idx, idx_a) = c;
                    A.coeffRef (idx, idx_b) = c;


                    if (dim == DIM_3D)
                    {
                        idx_a = mesh->GetGlobalIndexOfPoint (i, j, k+1);
                        idx_b = mesh->GetGlobalIndexOfPoint (i, j, k-1);

                        A.coeffRef (idx, idx_a) = d;
                        A.coeffRef (idx, idx_b) = d;

                    }
                }
            }


    std::cout << INDENT << "The matrix on the cartesian grid has been constructed." << std::endl;


    // On construit la partie supérieure de la matrice à partir de la liste de triplets


    // On lui ajoute sa transposée pour avoir la matrice complète

    //    A = A.selfadjointView<Eigen::Upper> ();
    //    A += A.transpose ();

    //    std::cout << A << std::endl;

    // On met à jour les interactions

    for (int k = 0; k < N-Ngrid; k++)
    {
            Point* P_k = mesh->GetPoint (k + Ngrid);
            std::vector <Point*> Neighbours = P_k->GetListNeighbours ();

            Point* P_l = Neighbours [0]; // Un voisin en direction gamma
            Point* P_r = Neighbours [1]; // L'autre voisin en direction gamma

            int l = P_l->GetGlobalIndex ();
            int r = P_r->GetGlobalIndex ();

            AXIS_LABEL gamma;
            Point Diff = *P_r - *P_l;

            if (Diff == Point (Diff.x, 0, 0)) // Les coordonnées x diffèrent, donc direction x
                gamma = AXIS_X;
            else if (Diff == Point (0, Diff.y, 0)) // Les coordonnées y diffèrent, donc direction y
                gamma = AXIS_Y;
            else // Les coordonnées z diffèrent, donc direction z
                gamma = AXIS_Z;


            A.coeffRef (l,r) = 0.;
            A.coeffRef (r,l) = 0.;

            Actualise_Ligne (A, P_k, gamma);
            Actualise_Ligne (A, P_l, gamma);
            Actualise_Ligne (A, P_r, gamma);


        std::cout << "\r" << INDENT << "Insertion of border point " << k+1 << "/" << N-Ngrid << "." << std::flush;
    }

    std::cout << "\r" << INDENT << "All border point are in the matrix now." << std::endl;
    A.pruned ();
    std::cout << INDENT << "Numeric zeros have been deleted." << std::endl;

    std::cout << std::endl;

    return A;
}


void Actualise_Ligne (Matrix &A, Point* P_m, AXIS_LABEL gamma)
{
    Sort_Neighbours (P_m);
    int m = P_m->GetGlobalIndex ();

    std::vector<Point*> Neighbours = P_m->GetListNeighbours ();

//    std::cout << "Point : " << m << ", neigh : " << Neighbours.size () << " AXE : " << gamma << std::endl;
    int beg = Remainder (int (2 * size_t(gamma)), int (Neighbours.size ()));
    int end = Remainder (int (2 * size_t(gamma) + 1), int (Neighbours.size ()));

    Point* P_l = Neighbours.at (size_t (beg));
    Point* P_r = Neighbours.at (size_t (end));

    int l = P_l->GetGlobalIndex (); // Voisin de "gauche" en direction gamma
    int r = P_r->GetGlobalIndex (); // Voisin de "droite" en direction gamma

    double dist_l = EuclidianDist (*P_m, *P_l);
    double dist_r = EuclidianDist (*P_m, *P_r);
    double moy = (dist_l + dist_r) / 2.;

    //    std::cout << "\tm : " << m << std::endl;
    //    std::cout << "\t\tL : " << l << "dist :" << dist_l << std::endl;
    //    std::cout << "\t\tR : " << r << "dist :" << dist_r<< std::endl;


    A.coeffRef (m,l) = 1. / (moy * dist_l);
    A.coeffRef (m,r) = 1. / (moy * dist_r);
    A.coeffRef (m,m) = 0.;
    A.coeffRef (m,m) = - A.row (m).sum ();
}


void Sort_Neighbours (Point* P)
{
    // On crée un vecteur NewList qui contiendra les voisins triés
    std::vector<Point*> Neighbours = P->GetListNeighbours ();
    size_t size = Neighbours.size ();
    std::vector<Point*> NewList (size);

    // Si size = 2 : rien à faire, les voisins sont déjà regroupés en direction x

    // Si size = 4 : on met d'abord ceux en direction x, puis y
    if (size == 4) {

        size_t i = 0;
        size_t j = 0;

        for (Point* V : Neighbours) {

            Point Diff = *V - *P;

            if (Diff == Point (Diff.x, 0, 0)) { // V a même y que P, donc direction x
                NewList.at (i) = V;
                i++;
            }
            else if (Diff == Point (0, Diff.y, 0)) { // // V a même x que P, donc direction y
                NewList.at (2 + j) = V;
                j++;
            }
            else
            {
                std::cout << INDENT << "ERROR::Sort_Neighbours it seems that the points are not aligned (point : " << P->GetGlobalIndex () << ")." << std::endl;
            }
        }
    }
    else if (size == 6) // Si size = 6 : on met d'abord ceux en direction x, puis y, puis z
    {
        size_t i = 0;
        size_t j = 0;
        size_t k = 0;

        for (Point* V : Neighbours) {

            Point Diff = *V - *P;

            if (Diff == Point (Diff.x, 0, 0)) { // V a même y et z que P, donc direction x
                NewList.at (i) = V;
                i++;
            }
            else if (Diff == Point (0, Diff.y, 0)) { // // V a même x et z que P, donc direction y
                NewList.at (2+j) = V;
                j++;
            }
            else if (Diff == Point (0, 0, Diff.z)) { // // V a même x et y que P, donc direction z
                NewList.at (4 + k) = V;
                k++;
            } else
            {
                std::cout << INDENT << "ERROR::Sort_Neighbours it seems that the points are not aligned (point : " << P->GetGlobalIndex () << ")." << std::endl;
            }
        }
    }

    P->SetListNeighbours (NewList);
}


void InsertBeta (Mesh* mesh, Matrix* A, Vector* beta_vec)
{
    int N = mesh->GetNumberOfTotalPoints ();

    size_t Ns = size_t(N);
    std::vector<double*> diags (Ns);

    for (int i = 0; i < N; ++i)
    {
        diags.at (size_t (i)) = &A->coeffRef (i, i);
        *diags.at (size_t (i)) = 0.;
    }

    for (int idx = 0; idx < N; ++idx)
    {
        Point* p = mesh->GetPoint (idx);

        double beta_p = beta_vec->coeffRef (idx);

        auto neighbours = p->GetListNeighbours ();

        double* diag = diags.at (size_t (idx));

        for (Point* p_neigh : neighbours)
        {
            int idx_neigh = p_neigh->GetGlobalIndex ();

            double beta_neigh = beta_vec->coeffRef (idx_neigh);

            double* coeff = &A->coeffRef (idx, idx_neigh);

            *coeff *= (beta_p + beta_neigh) / 2.;
            *diag -= *coeff;
        }
    }

    return;
}
