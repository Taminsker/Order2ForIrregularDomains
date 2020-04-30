#include "wfield.h"



Field GetWField(Mesh* mesh, Vector* phi, Vector* sol, std::vector<int>* idxsBorder, double h0)
{
    Field field = BuildWOnBorder (mesh, phi, sol, idxsBorder, h0);
    ExtendWToAllDomain (&field, mesh,  idxsBorder);

    return field;
}

Field BuildWOnBorder (Mesh* mesh, Vector* phi, Vector* sol, std::vector<int>* idxsBorder, double h0)
{
    std::cout << "# Build W vectors on border." << std::endl;

    size_t NumPoints = size_t (mesh->GetNumberOfTotalPoints ());

    DIM dim = mesh->GetDimension ();

    Field field;

//    std::cout << "W size : " << field.W.size () << " ikdsghabsh" << std::endl;
//    std::cout << field.W << std::endl;

//    field.W.resize (NumPoints);
//    field.Normals.resize (NumPoints);
//    field.GradPhi.resize (NumPoints);
//    field.GradTemperature.resize (NumPoints);

//    std::cout << "gradphi size : " << field.GradPhi.size ()<< std::endl;
//    std::cout << "W size : " << field.W.size ()<< std::endl;
//    std::cout << "Normals size : " << field.Normals.size ()<< std::endl;
//    std::cout << "GradTemperature size : " << field.GradTemperature.size ()<< std::endl;

    for (size_t i = 0; i < NumPoints; ++i)
    {
//        std::cout << "i : " << i << std::endl;
        field.W.push_back (new Point());
        field.Normals.push_back (new Point());
        field.GradPhi.push_back (new Point());
        field.GradTemperature.push_back (new Point());
    }

    for (int idx : *idxsBorder)
    {
        std::cout << "\r" << INDENT << "Point " << idx << std::flush;

        Point* p = mesh->GetPoint (idx);

        std::vector<Point*> listNeigh_p = p->GetListNeighbours ();

        size_t NumNeigh = listNeigh_p.size ();
        std::vector<double> Quantity (NumNeigh, 0.);

        for (size_t i = 0; i < NumNeigh; ++i)
        {
            Point* p_n = listNeigh_p.at (i);

            Point dist;

            std::vector<Point*> listNeigh_p_n = p_n->GetListNeighbours ();

            // Gradient de Phi en p_n
            Point* GradPhi = field.GradPhi.at (size_t (p_n->GetGlobalIndex ()));
            // Gradient de T en p_n
            Point* GradT = field.GradTemperature.at (size_t (p_n->GetGlobalIndex ()));
            Point* Normal = field.Normals.at(size_t (p_n->GetGlobalIndex ()));
            GradPhi->operator= ({0., 0., 0.});
            GradT->operator= ({0., 0., 0.});
            Normal->operator= ({0., 0., 0.});


            for (Point* v : listNeigh_p_n)
            {
                Point diff = *p_n - *v;

                double signe = (*p <= *v ? 1: -1);

                if (diff == Point(diff.x, 0, 0))
                {

                    GradPhi->x = GradPhi->x + signe * double(phi->coeff(v->GetGlobalIndex ()));
                    GradT->x = GradT->x + signe * double(sol->coeff(v->GetGlobalIndex ()));
                    dist.x = dist.x + fabs(diff.x);

                } else if (diff == Point(0, diff.y, 0) && dim >= DIM_2D)
                {
                    GradPhi->y = GradPhi->y + signe * double(phi->coeff(v->GetGlobalIndex ()));
                    GradT->y = GradT->y + signe * double(sol->coeff(v->GetGlobalIndex ()));
                    dist.y = dist.y + fabs(diff.y);

                } else if (diff == Point(0, 0, diff.z) && dim == DIM_3D)
                {
                    GradPhi->z = GradPhi->z + signe * (phi->coeff(v->GetGlobalIndex ()));
                    GradT->z = GradT->y + signe * (sol->coeff(v->GetGlobalIndex ()));
                    dist.z = dist.z + fabs(diff.z);
                }
                else
                {
                    std::cerr << INDENT << "ERROR Build on border :" << std::endl;
                    std::cerr << INDENT << "Point (" << p_n->GetGlobalIndex () << ")\t-\t" << *p_n << std::endl;
                    std::cerr << INDENT << "and Point (" << v->GetGlobalIndex () << ")\t-\t" << *v << std::endl;
                    std::cerr << INDENT << "Are not aligned but neighbours..." << std::endl;
                    exit(0);
                }
            }

            GradPhi->x = GradPhi->x / dist.x;
            if (dim >= DIM_2D)
                GradPhi->y = GradPhi->y / dist.y;
            if (dim == DIM_3D)
                GradPhi->z = GradPhi->z / dist.z;

            *Normal = *GradPhi / std::sqrt (*GradPhi|*GradPhi);

            GradT->x = GradT->x / dist.x;
            if (dim >= DIM_2D)
                GradT->y = GradT->y / dist.y;
            if (dim == DIM_3D)
                GradT->z = GradT->z / dist.z;

            Quantity.at (i) = (*GradT|*Normal);
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

        Point* normal_min = field.Normals.at (size_t (p_min->GetGlobalIndex ()));
        Point* normal_max = field.Normals.at (size_t (p_max->GetGlobalIndex ()));

        Point Normal_p = dist_min * *normal_min + dist_max * *normal_max;
        Normal_p = Normal_p / (dist_min + dist_max);
        *field.Normals.at (size_t (p->GetGlobalIndex ())) = Normal_p;

        *field.W.at (size_t (idx)) = S * Normal_p;
    }

    std::cout << "\r" << INDENT << "Build W on border is done.       " << std::endl;

    return field;
}


void ExtendWToAllDomain (Field* field, Mesh* mesh, std::vector<int>* idxsBorder)
{
    std::cout << "# Extend W to all domain." << std::endl;


    std::vector<int> indexes = mesh->GetListOfIndexPoints ();

    // LapW = 0 avec W = W_bordIdxs

    DIM dim = mesh->GetDimension ();
    int NumPoints = mesh->GetNumberOfTotalPoints ();

    Matrix A = Laplacian (mesh);
    RemovePeriodicity (mesh, &A);

    Vector bx(Vector::Zero (NumPoints));
    Vector by = bx;
    Vector bz = bx;

    for (int i : *idxsBorder)
        A.row (i) *= 0.;

    A = A.transpose ();

    for (int i : *idxsBorder)
    {
        // On déplace au 2nd membre les apparitions de P_i avec valeur imposée g(P_i)

        Point* w = field->W.at (size_t(i));

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

    ImposeZeroDirichletExtBorder (mesh, &A, &bx, &by, &bz);

    Wx = Solve (A, bx, IMPLICIT);

    if (dim >= DIM_2D)
        Wy = Solve (A, by, IMPLICIT);

    if (dim == DIM_3D)
        Wz = Solve (A, bz, IMPLICIT);


    size_t NumCartesian = size_t (mesh->GetNumberOfCartesianPoints ());

    for (size_t i = 0; i < NumCartesian; ++i)
    {
        delete field->W.at (i);

        Point* w = new Point();

        w->x = Wx.coeff(int(i));
        if (dim >= DIM_2D)
            w->y = Wy.coeff (int(i));
        if (dim == DIM_3D)
            w->z = Wz.coeff (int(i));

        field->W.at (i) = w;
    }

    std::cout << INDENT << "Extend W is done.       " << std::endl;

    return;
}
