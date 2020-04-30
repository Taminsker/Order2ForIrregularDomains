#include "impose.h"

void ImposeDirichlet (Mesh * mesh,
                      Matrix * A,
                      Vector* secondMember,
                      double (*g) (Point, double),
                      std::vector <int> listIndex,
                      double t)
{
    for (int i : listIndex)
    {
        A->row (i) *= 0.;
    } // On met la i-ème ligne à 0

    *A = A->transpose (); // Pour pouvoir manipuler les colonnes de A

    for (int i : listIndex) {
        // On déplace au 2nd membre les apparitions de P_i avec valeur imposée g(P_i)

        *secondMember -= g(*mesh->GetPoint (i), t) * A->row (i).transpose ();

        A->row (i) *= 0.;
        A->coeffRef (i,i) = 1.;

        secondMember->coeffRef (i) = g(*mesh->GetPoint (i), t);
    }

    *A = A->transpose ().pruned ();

    //    std::cout << *A << std::endl;

    return;
}


void ImposeDirichlet (Mesh * mesh,
                      Matrix * A,
                      Vector* secondMember,
                      Vector* g_list, // construit avec Funtovec 1
                      std::vector <int> listIndex,
                      double t)
{
    (void) mesh;
    (void)t;

    for (int i : listIndex) {A->row (i) *= 0.;} // On met la i-ème ligne à 0

    *A = A->transpose (); // Pour pouvoir manipuler les colonnes de A

    for (int i : listIndex) {
        // On déplace au 2nd membre les apparitions de P_i avec valeur imposée g(P_i)

        *secondMember -= g_list->coeff (i) * A->row (i).transpose ();

        A->row (i) *= 0.;
        A->coeffRef (i,i) = 1.;

        secondMember->coeffRef (i) = g_list->coeff (i);
    }

    *A = A->transpose ().pruned ();

    //    std::cout << *A << std::endl;

    return;
}

void ImposeNeumann (  Mesh *mesh,
                      Matrix * A,
                      Vector* scMember,
                      double (*g) (Point, double),
                      Point (*phigrad) (Point, double),
                      std::vector <int> listOfIndexes,
                      INTERPOLATION_NORMAL interpolationType,
                      double t)
{
    int numPoints = int(listOfIndexes.size ());

    std::cout << "# ImposeNeumann implicit function." << std::endl;
    std::cout << INDENT << "There are " << numPoints << " points selected." << std::endl;
    std::cout << INDENT << "The expected degree of interpolation is " << interpolationType << std::endl;

    // Créer le vecteur de conditions aux bords
    Vector cond_border (mesh->GetNumberOfTotalPoints ());

    DIM dim = mesh->GetDimension ();

    for (int idx : listOfIndexes)
    {
        Point* p = mesh->GetPoint (idx);
        auto neigh = p->GetListNeighbours ();

        size_t numOfNeigh = neigh.size ();

        if (numOfNeigh != 2 * size_t (dim) && idx < mesh->GetNumberOfCartesianPoints ())
        {
            std::cout << INDENT << "ERROR::ImposeNeumann the point " << idx << " has " << numOfNeigh << std::flush;
            std::cout << "\t\tnumber of neighbors expected : " << 2 << " or " << 2 * dim << "." << std::endl;

            continue;
        }

        std::vector<int> DF_index = {};
        std::vector<double> DF_coeff = {};

        switch (interpolationType) {
        case DEGRE_1:
            DF_index = {0, 1};
            DF_coeff = {-1, 1};
            break;
        case DEGRE_2:
            DF_index = {-1, 1};
            DF_coeff = {-1./2., 1./2.};
            break;
        case DEGRE_3:
            DF_index = {0, 1, 2, 3};
            DF_coeff = {-11./6., 3., -3./2., 1./3};
            break;
        case DEGRE_4:
            DF_index = {-2, -1, 1, 2};
            DF_coeff = {1./12., -2./3., 2./3., -1./12.};
            break;
        }

        Point normal = phigrad (*p, t);

        bool NotOnTheGrid = (numOfNeigh != 2 * size_t (dim));

        std::vector<Point*> vec_point = {};
        std::vector<double> vec_coeff = {};
        std::vector<double> h = {normal.x * mesh->Get_hx (),
                                 normal.y * mesh->Get_hy (),
                                 normal.z * mesh->Get_hz ()};

        Point* neigh1 = neigh.at (0);
        Point* neigh2 = neigh.at (1);

        double rapport_1 = EuclidianDist (*p, *neigh1) / EuclidianDist (*neigh1, *neigh2);
        double rapport_2 = EuclidianDist (*p, *neigh2) / EuclidianDist (*neigh1, *neigh2);

        LocalIndexes l1 = mesh->GetLocalIndexesOfPoint (neigh1->GetGlobalIndex ());
        LocalIndexes l2 = mesh->GetLocalIndexesOfPoint (neigh2->GetGlobalIndex ());
        LocalIndexes lp = mesh->GetLocalIndexesOfPoint (p->GetGlobalIndex ());

        for (size_t dimension = 0; dimension < size_t(dim); dimension++)
        {
            std::vector<int> a = {int(dimension == 0), int(dimension == 1), int(dimension == 2)};

            for (size_t looper = 0; looper < DF_index.size (); ++looper)
            {
                int i = a.at (0) * DF_index.at (looper);
                int j = a.at (1) * DF_index.at (looper);
                int k = a.at (2) * DF_index.at (looper);

                if (NotOnTheGrid)
                {
                    vec_point.push_back (mesh->GetPoint (l1.i + i, l1.j + j, l1.k + k));
                    vec_coeff.push_back (DF_coeff.at (looper) * h.at (dimension) * rapport_1);

                    vec_point.push_back (mesh->GetPoint (l2.i + i, l2.j + j, l2.k + k));
                    vec_coeff.push_back (DF_coeff.at (looper) * h.at (dimension) * rapport_2);
                }
                else
                {
                    vec_point.push_back (mesh->GetPoint (lp.i + i, l1.j + j, l1.k + k));
                    vec_coeff.push_back (DF_coeff.at (looper) * h.at (dimension));
                }
            }
        }

        A->row (idx) *= 0;

        scMember->coeffRef (idx) = g(*p, t);

        for (size_t i = 0; i < vec_coeff.size (); ++i)
        {
            int col = vec_point.at (i)->GetGlobalIndex ();
            A->coeffRef (idx, col) += vec_coeff.at (i);
        }

        std::cout << "\r" << INDENT << "The point " << idx << " was imposed.          " << std::flush;
    }

    std::cout << "\r" << INDENT << "All points have been imposed.                                 " << std::endl;

    A->pruned ();

    return;
}

void ImposeNeumann (  Mesh *mesh,
                      Matrix * A,
                      Vector* scMember,
                      double (*g) (Point, double),
                      Point (*phi) (double, double),
                      Point (*phigrad) (double, double),
                      std::vector <int> listOfIndexes,
                      INTERPOLATION_NORMAL interpolationType,
                      double t,
                      double a,
                      double b)
{

    Point centre;
    for (int i = 0; i < 300; ++i)
        centre = centre + phi(a + (b-a)/300., t);

    centre = centre / 300.;

    int numPoints = int(listOfIndexes.size ());

    std::cout << "# ImposeNeumann parametric function." << std::endl;
    std::cout << INDENT << "There are " << numPoints << " points selected." << std::endl;
    std::cout << INDENT << "The expected degree of interpolation is " << interpolationType << std::endl;

    // Créer le vecteur de conditions aux bords
    Vector cond_border (mesh->GetNumberOfTotalPoints ());

    DIM dim = mesh->GetDimension ();

    for (int idx : listOfIndexes)
    {
        Point* p = mesh->GetPoint (idx);
        auto neigh = p->GetListNeighbours ();

        size_t numOfNeigh = neigh.size ();

        if (numOfNeigh != 2 * size_t (dim) && idx < mesh->GetNumberOfCartesianPoints ())
        {
            std::cout << INDENT << "ERROR::ImposeNeumann the point " << idx << " has " << numOfNeigh << std::flush;
            std::cout << "\t\tnumber of neighbors expected : " << 2 << " or " << 2 * dim << "." << std::endl;

            continue;
        }

        std::vector<int> DF_index = {};
        std::vector<double> DF_coeff = {};

        switch (interpolationType) {
        case DEGRE_1:
            DF_index = {0, 1};
            DF_coeff = {-1, 1};
            break;
        case DEGRE_2:
            DF_index = {-1, 1};
            DF_coeff = {-1./2., 1./2.};
            break;
        case DEGRE_3:
            DF_index = {0, 1, 2, 3};
            DF_coeff = {-11./6., 3., -3./2., 1./3};
            break;
        case DEGRE_4:
            DF_index = {-2, -1, 1, 2};
            DF_coeff = {1./12., -2./3., 2./3., -1./12.};
            break;
        }

        Point normal = phigrad (std::atan ((p->y - centre.y) / fabs(p->x - centre.x)), t);
        double temp = normal.x;
        normal.x = -1 * normal.y;
        normal.y = temp;
        normal = normal / std::sqrt(normal | normal);

        bool NotOnTheGrid = (numOfNeigh != 2 * size_t (dim));

        std::vector<Point*> vec_point = {};
        std::vector<double> vec_coeff = {};
        std::vector<double> h = {normal.x * mesh->Get_hx (),
                                 normal.y * mesh->Get_hy (),
                                 normal.z * mesh->Get_hz ()};

        Point* neigh1 = neigh.at (0);
        Point* neigh2 = neigh.at (1);

        double rapport_1 = EuclidianDist (*p, *neigh1) / EuclidianDist (*neigh1, *neigh2);
        double rapport_2 = EuclidianDist (*p, *neigh2) / EuclidianDist (*neigh1, *neigh2);

        LocalIndexes l1 = mesh->GetLocalIndexesOfPoint (neigh1->GetGlobalIndex ());
        LocalIndexes l2 = mesh->GetLocalIndexesOfPoint (neigh2->GetGlobalIndex ());
        LocalIndexes lp = mesh->GetLocalIndexesOfPoint (p->GetGlobalIndex ());

        for (size_t dimension = 0; dimension < size_t(dim); dimension++)
        {
            std::vector<int> a = {int(dimension == 0), int(dimension == 1), int(dimension == 2)};

            for (size_t looper = 0; looper < DF_index.size (); ++looper)
            {
                int i = a.at (0) * DF_index.at (looper);
                int j = a.at (1) * DF_index.at (looper);
                int k = a.at (2) * DF_index.at (looper);

                if (NotOnTheGrid)
                {
                    vec_point.push_back (mesh->GetPoint (l1.i + i, l1.j + j, l1.k + k));
                    vec_coeff.push_back (DF_coeff.at (looper) * h.at (dimension) * rapport_1);

                    vec_point.push_back (mesh->GetPoint (l2.i + i, l2.j + j, l2.k + k));
                    vec_coeff.push_back (DF_coeff.at (looper) * h.at (dimension) * rapport_2);
                }
                else
                {
                    vec_point.push_back (mesh->GetPoint (lp.i + i, l1.j + j, l1.k + k));
                    vec_coeff.push_back (DF_coeff.at (looper) * h.at (dimension));
                }
            }
        }

        A->row (idx) *= 0;

        scMember->coeffRef (idx) = g(*p, t);

        for (size_t i = 0; i < vec_coeff.size (); ++i)
        {
            int col = vec_point.at (i)->GetGlobalIndex ();
            A->coeffRef (idx, col) += vec_coeff.at (i);
        }

        std::cout << "\r" << INDENT << "The point " << idx << " was imposed.          " << std::flush;
    }

    std::cout << "\r" << INDENT << "All points have been imposed.                                 " << std::endl;

    A->pruned ();

    return;
}

void ImposeZeroDirichletExtBorder (Mesh* mesh, Matrix* A, Vector* a1, Vector* a2, Vector* a3)
{

    std::cout << "# Impose 0 on exterieur border " << std::endl;

    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    //    int N = mesh->GetNumberOfCartesianPoints ();
    DIM dim = mesh->GetDimension ();

    std::vector<int> indexes = {};

    if (dim == DIM_1D)
    {
        indexes = {0, Nx-1};

    } else if (dim == DIM_2D)
    {
        for (int i = 0; i < Nx; ++i)
        {
            indexes.push_back (mesh->GetGlobalIndexOfPoint (i, 0));
            indexes.push_back (mesh->GetGlobalIndexOfPoint (i, Ny-1));
        }

        for (int j = 0; j < Ny; ++j)
        {
            // Interaction gauche droite

            indexes.push_back (mesh->GetGlobalIndexOfPoint (0, j));
            indexes.push_back (mesh->GetGlobalIndexOfPoint (Nx-1, j));
        }
    } else
    {

        for (int k = 0; k < Nz; ++k)
        {
            for (int i = 0; i < Nx; ++i)
            {
                indexes.push_back (mesh->GetGlobalIndexOfPoint (i, 0, k));
                indexes.push_back (mesh->GetGlobalIndexOfPoint (i, Ny-1, k));
            }

            for (int j = 0; j < Ny; ++j)
            {
                indexes.push_back (mesh->GetGlobalIndexOfPoint (0, j, k));
                indexes.push_back (mesh->GetGlobalIndexOfPoint (Nx-1, j, k));
            }
        }

        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                indexes.push_back (mesh->GetGlobalIndexOfPoint (i, j, 0));
                indexes.push_back (mesh->GetGlobalIndexOfPoint (i, j, Nz-1));
            }

            for (int k = 0; k < Nz; ++k)
            {
                indexes.push_back (mesh->GetGlobalIndexOfPoint (i, 0, k));
                indexes.push_back (mesh->GetGlobalIndexOfPoint (i, Ny-1, k));
            }
        }

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                indexes.push_back (mesh->GetGlobalIndexOfPoint (i, j, 0));
                indexes.push_back (mesh->GetGlobalIndexOfPoint (i, j, Nz-1));
            }

            for (int k = 0; k < Nz; ++k)
            {
                indexes.push_back (mesh->GetGlobalIndexOfPoint (0, j, k));
                indexes.push_back (mesh->GetGlobalIndexOfPoint (Nx-1, j, k));
            }
        }
    }

    Vector Zero;
    Zero.resize (mesh->GetNumberOfTotalPoints ());
    Zero.setZero ();
    auto C1 = *A;

    if (a1 != nullptr)
        ImposeDirichlet (mesh, &C1, a1, &Zero, indexes);

    if (a2 != nullptr)
    {
        C1 = *A;
        ImposeDirichlet (mesh, &C1, a2, &Zero, indexes);
    }


    if (a3 != nullptr)
    {
        C1 = *A;
        ImposeDirichlet (mesh, &C1, a3, &Zero, indexes);
    }

    *A = C1.pruned ();

    std::cout << std::endl;

    return;
}
