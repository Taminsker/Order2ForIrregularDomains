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

    // On définit un alias pour les triplets

    typedef Eigen::Triplet<double> Triplet;

    // On définit la liste de triplets et on prépare la place

    std::vector<Triplet> ListOfTriplets;
    ListOfTriplets.reserve(size_t (7 * N));

    // On remplit la liste de triplets

    //    for (int i = 0; i < Ngrid; ++i) {

    //        int surDiag1 = i + 1;
    //        int surDiag2 = i + Nx;
    //        int surDiag3 = i + (Nx * Ny);

    //        // La diagonale principale des a

    //        ListOfTriplets.push_back(Triplet(i,i,a));

    //        // La sur-diagonale des b, décalée de 1 (si possible)

    //        if (surDiag1 < Ngrid)
    //            ListOfTriplets.push_back(Triplet(i,surDiag1,b));

    //        // La sur-diagonale des c, décalée de Nx (si possible)

    //        if (surDiag2 < Ngrid)
    //            ListOfTriplets.push_back(Triplet(i,surDiag2,c));

    //        // La sur-diagonale des d, décalée de Nx * Ny (si possible)

    //        if (surDiag3 < Ngrid)
    //            ListOfTriplets.push_back(Triplet(i,surDiag3,d));

    //    }

    DIM dim = mesh->GetDimension ();

    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
            {

                int idx = mesh->GetGlobalIndexOfPoint (i, j, k);

                ListOfTriplets.push_back (Triplet (idx, idx, a));

                int idx_a = mesh->GetGlobalIndexOfPoint (i+1, j, k);
                int idx_b = mesh->GetGlobalIndexOfPoint (i-1, j, k);

                ListOfTriplets.push_back (Triplet (idx,idx_a, b));
                ListOfTriplets.push_back (Triplet (idx,idx_b, b));

                if (dim == DIM_2D || dim == DIM_3D)
                {
                    idx_a = mesh->GetGlobalIndexOfPoint (i, j+1, k);
                    idx_b = mesh->GetGlobalIndexOfPoint (i, j-1, k);

                    ListOfTriplets.push_back (Triplet (idx,idx_a, c));
                    ListOfTriplets.push_back (Triplet (idx,idx_b, c));

                    if (dim == DIM_3D)
                    {
                        idx_a = mesh->GetGlobalIndexOfPoint (i, j, k+1);
                        idx_b = mesh->GetGlobalIndexOfPoint (i, j, k-1);

                        ListOfTriplets.push_back (Triplet (idx,idx_a, d));
                        ListOfTriplets.push_back (Triplet (idx,idx_b, d));
                    }
                }
            }


    std::cout << INDENT << "The matrix on the cartesian grid has been constructed." << std::endl;


    // On construit la partie supérieure de la matrice à partir de la liste de triplets

    Matrix A (N, N);
    A.setFromTriplets(ListOfTriplets.begin(), ListOfTriplets.end());

    // On lui ajoute sa transposée pour avoir la matrice complète

    //    A = A.selfadjointView<Eigen::Upper> ();
    //    A += A.transpose ();

//    std::cout << A << std::endl;

    // On met à jour les interactions

    for (int k = 0; k < N; k++)
    {
        if (k >= Ngrid)
        {
            Point* P_k = mesh->GetPoint (k);
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
        }

        std::cout << "\r" << INDENT << "Insertion of border point " << k+1 << "/" << N << "." << std::flush;
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

    Point* P_l = Neighbours.at (2 * size_t(gamma));
    Point* P_r = Neighbours.at (2 * size_t(gamma) + 1);

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
