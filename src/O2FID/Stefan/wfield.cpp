#include "wfield.h"



std::vector<Point*> GetWField(Mesh* mesh, Vector* phi, Vector* sol, std::vector<int>* idxsBorder, double h0)
{
    auto W = BuildWOnBorder (mesh, phi, sol, idxsBorder, h0);
    ExtendWToAllDomain (mesh, &W, idxsBorder);
    return W;
}

std::vector<Point*> BuildWOnBorder (Mesh* mesh, Vector* phi, Vector* sol, std::vector<int>* idxsBorder, double h0)
{
    std::cout << "# Build W vectors on border." << std::endl;

    size_t NumPoints = size_t (mesh->GetNumberOfTotalPoints ());

    DIM dim = mesh->GetDimension ();

    std::vector<Point*> W (NumPoints);

    for (int idx : *idxsBorder)
    {
        std::cout << "\r" << INDENT << "Point " << idx << std::flush;

        Point* p = mesh->GetPoint (idx);

        std::vector<Point*> listNeigh_p = p->GetListNeighbours ();

        size_t NumNeigh = listNeigh_p.size ();
        std::vector<double> Quantity (NumNeigh, 0.);
        std::vector<Point> Normals (NumNeigh, Point());

        for (size_t i = 0; i < NumNeigh; ++i)
        {
            Point* p_n = listNeigh_p.at (i);

            Point dist;

            std::vector<Point*> listNeigh_p_n = p_n->GetListNeighbours ();

            // Gradient de Phi en p_n
            Point GradPhi;
            // Gradient de T en p_n
            Point GradT;

            for (Point* v : listNeigh_p_n)
            {
                Point diff = *p_n - *v;

                double signe = (*p <= *v ? 1: -1);

                if (diff == Point(diff.x, 0, 0))
                {
                    GradPhi.x += signe * phi->coeff(v->GetGlobalIndex ());
                    GradT.x += signe * sol->coeff(v->GetGlobalIndex ());

                    dist.x += fabs(diff.x);

                } else if (diff == Point(0, diff.y, 0) && dim >= DIM_2D)
                {
                    GradPhi.y += signe * phi->coeff(v->GetGlobalIndex ());
                    GradT.y += signe * sol->coeff(v->GetGlobalIndex ());

                    dist.y += fabs(diff.y);

                } else if (diff == Point(0, 0, diff.z) && dim == DIM_3D)
                {
                    GradPhi.z += signe * phi->coeff(v->GetGlobalIndex ());
                    GradT.z += signe * sol->coeff(v->GetGlobalIndex ());

                    dist.z += fabs(diff.z);
                }
                else
                {
                    std::cerr << "ERROR" << std::endl;
                    exit(0);
                }
            }

            GradPhi.x = GradPhi.x / dist.x;
            if (dim >= DIM_2D)
                GradPhi.y = GradPhi.y / dist.y;
            if (dim == DIM_3D)
                GradPhi.z = GradPhi.z / dist.z;

            Normals.at (i) = GradPhi / std::sqrt (GradPhi|GradPhi);

            GradT.x = GradT.x / dist.x;
            if (dim >= DIM_2D)
                GradT.y = GradT.y / dist.y;
            if (dim == DIM_3D)
                GradT.z = GradT.z / dist.z;

            Quantity.at (i) = (GradT|Normals.at (i));
        }

        size_t idx_min = 0;
        size_t idx_max = idx_min;

        double phi_min = sol->coeff (listNeigh_p.front ()->GetGlobalIndex ());
        double phi_max = phi_min;

        for (size_t i = 0; i < NumNeigh; ++i)
        {
            Point* p_n = listNeigh_p.at (i);

            double value = sol->coeff (p_n->GetGlobalIndex ());
            if (value < phi_min)
            {
                phi_min = value;
                idx_min = i;
            }

            if (value > phi_max)
            {
                phi_max = value;
                idx_max = i;
            }
        }

        double S = -1. / h0 * (Quantity.at (idx_min) - Quantity.at (idx_max));

        Point* p_min = listNeigh_p.at (idx_min);
        Point* p_max = listNeigh_p.at (idx_max);

        double dist_min = EuclidianDist (*p, *p_min);
        double dist_max = EuclidianDist (*p, *p_max);

        Point normal_min = Normals.at (idx_min);
        Point normal_max = Normals.at (idx_max);

        Point Normal_p = dist_min * normal_min + dist_max * normal_max;
        Normal_p = Normal_p / (dist_min + dist_max);

        W.at (size_t (idx)) = new Point(S * Normal_p);
    }

    std::cout << "\r" << INDENT << "Build W on border is done.       " << std::endl;

    return W;
}


void ExtendWToAllDomain (Mesh* mesh, std::vector<Point*>* W, std::vector<int>* idxsBorder)
{
    std::vector<int> indexes = mesh->GetListOfIndexPoints ();

    // LapW = 0 avec W = W_bordIdxs

    DIM dim = mesh->GetDimension ();
    int NumPoints = mesh->GetNumberOfTotalPoints ();

    Matrix A = Laplacian (mesh);

    Vector bx(Vector::Zero (NumPoints));
    Vector by = bx;
    Vector bz = bx;

    for (int i : *idxsBorder)
        A.row (i) *= 0.;

    A = A.transpose ();

    for (int i : *idxsBorder)
    {
        // On déplace au 2nd membre les apparitions de P_i avec valeur imposée g(P_i)

        Point* w = W->at (size_t(i));

        bx -= w->x * A.row (i).transpose ();

        if (dim >= DIM_2D)
            by -= w->y * A.row (i).transpose ();
        if (dim == DIM_3D)
            bz -= w->z * A.row (i).transpose ();

        A.row (i) *= 0.;
        A.coeffRef (i,i) = 1.;

        bx.coeffRef (i) = w->x;
        if (dim >= DIM_2D)
            by.coeffRef (i) = w->y;
        if (dim == DIM_3D)
            bz.coeffRef (i) = w->z;
    }

    A = A.transpose ().pruned ();

    Vector Wx, Wy, Wz = bx;

    RemovePeriodicity (mesh, &A);

    Wx = Solve (A, bx, IMPLICIT);

    if (dim >= DIM_2D)
        Wy = Solve (A, by, IMPLICIT);

    if (dim == DIM_3D)
        Wz = Solve (A, bz, IMPLICIT);

//    AutoClearVector (W);

    size_t NumCartesian = size_t (mesh->GetNumberOfCartesianPoints ());

//    *W = std::vector<Point*>(NumCartesian);

    for (size_t i = 0; i < NumCartesian; ++i)
    {
        delete W->at (i);

        Point* w = new Point();

        w->x = Wx.coeff(int(i));
        if (dim >= DIM_2D)
            w->y = Wy.coeff (int(i));
        if (dim == DIM_3D)
            w->z = Wz.coeff (int(i));

        W->at (i) = w;
    }

    std::cout << INDENT << "Extend W is done.       " << std::endl;

    return;
}
