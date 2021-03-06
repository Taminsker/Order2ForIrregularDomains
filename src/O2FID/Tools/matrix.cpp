/** @file matrix.cpp */


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

    if (hx > 1e-15) {b = 1. / (hx * hx);}
    if (hy > 1e-15) {c = 1. / (hy * hy);}
    if (hz > 1e-15) {d = 1. / (hz * hz);}
    a = - 2. * b - 2. * c - 2. * d;

    DIM dim = mesh->GetDimension ();

    Matrix A (N, N);
    A.reserve (Eigen::VectorXi::Constant(N, 15));

    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
            {
                //                std::cout << "i = " << i << std::endl;

                int idx = mesh->GetGlobalIndexOfPoint (i, j, k);

                int idx_L = mesh->GetGlobalIndexOfPoint (i-1, j, k); // Left
                int idx_R = mesh->GetGlobalIndexOfPoint (i+1, j, k); // Right

                A.coeffRef (idx,idx_L) = b;
                A.coeffRef (idx,idx_R) = b;

                if (dim == DIM_2D || dim == DIM_3D)
                {
                    int idx_F = mesh->GetGlobalIndexOfPoint (i, j-1, k); // Front
                    int idx_B = mesh->GetGlobalIndexOfPoint (i, j+1, k); // Back

                    // Si le point n'est pas sur un bord en direction y : Laplacien classique
                    A.coeffRef (idx,idx_F) = c;
                    A.coeffRef (idx,idx_B) = c;

                    if (dim == DIM_3D)
                    {
                        int idx_D = mesh->GetGlobalIndexOfPoint (i, j, k-1); // Down
                        int idx_U = mesh->GetGlobalIndexOfPoint (i, j, k+1); // Up

                        A.coeffRef (idx,idx_D) = d;
                        A.coeffRef (idx,idx_U) = d;

                    } // fin if dim 3
                } // fin if dim 2 ou 3

                A.coeffRef (idx,idx) = a;

            } // fin for i,j,k


    std::cout << INDENT << "The matrix on the cartesian grid has been constructed." << std::endl;


    // On construit la partie supérieure de la matrice à partir de la liste de triplets


    // On lui ajoute sa transposée pour avoir la matrice complète

    //    A = A.selfadjointView<Eigen::Upper> ();
    //    A += A.transpose ();

    //        std::cout << A << std::endl;

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

    //    moy = std::max (dist_l, dist_r); // pour rendre symétrique

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

                std::cout << INDENT << "diff = " << Diff << std::endl;
                std::cout << INDENT << "p = " << *P << std::endl;
                std::cout << INDENT << "p_n = " << *V << std::endl;

                exit(0);
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

void RemovePeriodicity (Mesh* mesh, Matrix* A)
{

    std::cout << "# Remove Periodicity..." << std::endl;

    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    //    int N = mesh->GetNumberOfCartesianPoints ();
    DIM dim = mesh->GetDimension ();


    if (dim == DIM_1D)
    {
        A->coeffRef (0, Nx-1) = 0.;
        A->coeffRef (Nx-1, 0) = 0.;

    } else if (dim == DIM_2D)
    {
        for (int i = 0; i < Nx; ++i)
        {
            // Interaction haut bas
            int idx1 = mesh->GetGlobalIndexOfPoint (i, 0);
            int idx2 = mesh->GetGlobalIndexOfPoint (i, Ny-1);

            A->coeffRef (idx1, idx2) = 0.;
            A->coeffRef (idx2, idx1) = 0.;
        }

        for (int j = 0; j < Ny; ++j)
        {
            // Interaction gauche droite

            int idx1 = mesh->GetGlobalIndexOfPoint (0, j);
            int idx2 = mesh->GetGlobalIndexOfPoint (Nx-1, j);

            A->coeffRef (idx1, idx2) = 0.;
            A->coeffRef (idx2, idx1) = 0.;
        }
    } else
    {
        int idx1;
        int idx2;

        for (int k = 0; k < Nz; ++k)
        {
            for (int i = 0; i < Nx; ++i)
            {
                idx1 = mesh->GetGlobalIndexOfPoint (i, 0, k);
                idx2 = mesh->GetGlobalIndexOfPoint (i, Ny-1, k);

                A->coeffRef (idx1, idx2) = 0.;
                A->coeffRef (idx2, idx1) = 0.;
            }

            for (int j = 0; j < Ny; ++j)
            {

                idx1 = mesh->GetGlobalIndexOfPoint (0, j, k);
                idx2 = mesh->GetGlobalIndexOfPoint (Nx-1, j, k);

                A->coeffRef (idx1, idx2) = 0.;
                A->coeffRef (idx2, idx1) = 0.;
            }
        }

        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                idx1 = mesh->GetGlobalIndexOfPoint (i, j, 0);
                idx2 = mesh->GetGlobalIndexOfPoint (i, j, Nz-1);

                A->coeffRef (idx1, idx2) = 0.;
                A->coeffRef (idx2, idx1) = 0.;
            }

            for (int k = 0; k < Nz; ++k)
            {
                idx1 = mesh->GetGlobalIndexOfPoint (i, 0, k);
                idx2 = mesh->GetGlobalIndexOfPoint (i, Ny-1, k);

                A->coeffRef (idx1, idx2) = 0.;
                A->coeffRef (idx2, idx1) = 0.;
            }
        }

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                idx1 = mesh->GetGlobalIndexOfPoint (i, j, 0);
                idx2 = mesh->GetGlobalIndexOfPoint (i, j, Nz-1);

                A->coeffRef (idx1, idx2) = 0.;
                A->coeffRef (idx2, idx1) = 0.;
            }

            for (int k = 0; k < Nz; ++k)
            {
                idx1 = mesh->GetGlobalIndexOfPoint (0, j, k);
                idx2 = mesh->GetGlobalIndexOfPoint (Nx-1, j, k);

                A->coeffRef (idx1, idx2) = 0.;
                A->coeffRef (idx2, idx1) = 0.;
            }
        }
    }

    *A = A->pruned ();

    std::cout << std::endl;

    return;
}

Matrix Gradient (Mesh* mesh, ORDERS order)
{
    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hy ();
    double hz = mesh->Get_hz ();

    double maxima = std::max (hx, std::max (hy, hz)) + 1e-15;

    int N = mesh->GetNumberOfTotalPoints ();
    DIM dim = mesh->GetDimension ();

    Matrix A(N, N);

    auto DF = DFStruct();

    std::vector<int> idxDf;
    std::vector<double> coeffDf;

    DFOrderBuild (1, order, &idxDf, &coeffDf);

    size_t SizeDf = idxDf.size ();

    A.reserve (Eigen::VectorXi::Constant(N, int(3 * idxDf.size ())));

    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                Point* p = mesh->GetPoint (i, j, k);
                int idxGlobal = p->GetGlobalIndex ();

                for (DIM d : {DIM_1D, DIM_2D, DIM_3D})
                {
                    if (d > dim)
                        break;

                    for (size_t s = 0; s < SizeDf; ++s)
                    {
                        double coeff = coeffDf.at (s);
                        int emp = idxDf.at (s);

                        int incI = int(d == DIM_1D);
                        int incJ = int(d == DIM_2D);
                        int incK = int(d == DIM_3D);

                        int I = i + incI * emp;
                        int J = j + incJ * emp;
                        int K = k + incK * emp;

                        Point* n = mesh->GetPoint (I, J, K);
                        int idx = n->GetGlobalIndex ();


                        switch (d)
                        {
                        case DIM_1D:
                            coeff /= hx;
                            break;
                        case DIM_2D:
                            coeff /= hy;
                            break;
                        case DIM_3D:
                            coeff /= hz;
                            break;
                        }

                        if (EuclidianDist (*p, *n) < maxima)
                            A.coeffRef (idxGlobal, idx) += coeff;
                    }
                }
            }
        }
    }

    auto listPointIdx = mesh->GetListOfIndexPoints ();

    int G = mesh->GetNumberOfCartesianPoints ();

    // Insertion des points de bords

    for (int idx : listPointIdx)
    {
        if (idx < G)
            continue;

        std::cout << "idx : " << idx << std::endl;

        Point* p = mesh->GetPoint (idx);
        auto neigh = p->GetListNeighbours ();

        Point* n1 = neigh.front ();
        Point* n2 = neigh.back ();
        Point* temp = n1;

        if (*n1 >= *n2)
        {
            n1 = n2;
            n2 = temp;
        }

        LocalIndexes loc1 = mesh->GetLocalIndexesOfPoint (n1->GetGlobalIndex ());
        LocalIndexes loc2 = mesh->GetLocalIndexesOfPoint (n2->GetGlobalIndex ());

        int Idiff = loc1.i - loc2.i; // Axe i
        int Jdiff = loc1.j - loc2.j; // Axe j
        int Kdiff = loc1.k - loc2.k; // Axe k

        if (Idiff != 0 && Jdiff != 0 && Kdiff != 0)
        {
            std::cout << INDENT << "Gradient detects two neighbours of point not aligned..." << std::endl;

            exit(0);
        }

        // ie n2 + {Idff, Jdiff, Kdiff} = n1
        // ie n1 - {Idff, Jdiff, Kdiff} = n2
        // et p est entre n1 et n2

        // Point n1
        int idxGlobal = n1->GetGlobalIndex ();
        int i = loc1.i;
        int j = loc1.j;
        int k = loc1.k;

        A.row (idxGlobal) *= 0.;

        for (size_t s = 0; s < SizeDf; ++s)
        {
            double coeff = coeffDf.at (s);
            int emp = idxDf.at (s);

            if (emp > Idiff || Idiff == 0)
            {
                int idxG = mesh->GetGlobalIndexOfPoint (i + emp, j, k);
                A.coeffRef (idxGlobal, idxG) = coeff / hx;
            } else if (emp < Idiff)
            {
                emp++;
                int idxG = mesh->GetGlobalIndexOfPoint (i + emp, j, k);
                A.coeffRef (idxGlobal, idxG) = coeff / hx;

            } else
            {
                A.coeffRef (idxGlobal, idx) = coeff / EuclidianDist (*n1, *p);
            }

            if (dim == DIM_2D || dim == DIM_3D)
            {
                if (emp > Jdiff || Jdiff == 0)
                {
                    int idxG = mesh->GetGlobalIndexOfPoint (i, j + emp, k);
                    A.coeffRef (idxGlobal, idxG) = coeff / hy;
                } else if (emp < Jdiff)
                {
                    emp++;
                    int idxG = mesh->GetGlobalIndexOfPoint (i, j + emp, k);
                    A.coeffRef (idxGlobal, idxG) = coeff / hy;

                } else
                {
                    A.coeffRef (idxGlobal, idx) = coeff / EuclidianDist (*n1, *p);
                }

                if (dim == DIM_3D)
                {
                    if (emp > Kdiff || Kdiff == 0)
                    {
                        int idxG = mesh->GetGlobalIndexOfPoint (i, j, k + emp);
                        A.coeffRef (idxGlobal, idxG) = coeff / hz;
                    } else if (emp < Kdiff)
                    {
                        emp++;
                        int idxG = mesh->GetGlobalIndexOfPoint (i, j, k + emp);
                        A.coeffRef (idxGlobal, idxG) = coeff / hz;

                    } else
                    {
                        A.coeffRef (idxGlobal, idx) = coeff / EuclidianDist (*n1, *p);
                    }
                }
            }
        }


        // Point n2
        idxGlobal = n2->GetGlobalIndex ();
        i = loc2.i;
        j = loc2.j;
        k = loc2.k;

        A.row (idxGlobal) *= 0.;

        for (size_t s = 0; s < SizeDf; ++s)
        {
            double coeff = coeffDf.at (s);
            int emp = idxDf.at (s);

            if (emp < -1 * Idiff || Idiff == 0)
            {
                int idxG = mesh->GetGlobalIndexOfPoint (i + emp, j, k);
                A.coeffRef (idxGlobal, idxG) = coeff / hx;
            } else if (emp > -1 * Idiff)
            {
                emp--;
                int idxG = mesh->GetGlobalIndexOfPoint (i + emp, j, k);
                A.coeffRef (idxGlobal, idxG) = coeff / hx;

            } else
            {
                A.coeffRef (idxGlobal, idx) = coeff / EuclidianDist (*n1, *p);
            }

            if (dim == DIM_2D || dim == DIM_3D)
            {
                if (emp < -1 * Jdiff || Jdiff == 0)
                {
                    int idxG = mesh->GetGlobalIndexOfPoint (i, j + emp, k);
                    A.coeffRef (idxGlobal, idxG) = coeff / hy;
                } else if (emp > -1 * Jdiff)
                {
                    emp--;
                    int idxG = mesh->GetGlobalIndexOfPoint (i, j + emp, k);
                    A.coeffRef (idxGlobal, idxG) = coeff / hy;

                } else
                {
                    A.coeffRef (idxGlobal, idx) = coeff / EuclidianDist (*n1, *p);
                }

                if (dim == DIM_3D)
                {
                    if (emp < -1 * Kdiff || Kdiff == 0)
                    {
                        int idxG = mesh->GetGlobalIndexOfPoint (i, j, k + emp);
                        A.coeffRef (idxGlobal, idxG) = coeff / hz;
                    } else if (emp > -1 * Kdiff)
                    {
                        emp--;
                        int idxG = mesh->GetGlobalIndexOfPoint (i, j, k + emp);
                        A.coeffRef (idxGlobal, idxG) = coeff / hz;

                    } else
                    {
                        A.coeffRef (idxGlobal, idx) = coeff / EuclidianDist (*n1, *p);
                    }
                }
            }
        }

        A.coeffRef (idx, idx) = 1.;
    }

    return A.pruned ();
}
