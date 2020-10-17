/** @file wfield.cpp*/

#include "wfield.h"
#include "interface.h"



Field* GetWField(Mesh* mesh, Vector* phi, Vector* sol, std::vector<int>* idxsBorder, double h0)
{
    Field* field = new Field();
    field->Normals = ComputeGradient (mesh, phi, true, ORDER_2_CENTRAL);
    field->GradPhi = ComputeGradient (mesh, phi, false, ORDER_2_CENTRAL);
    field->GradTemperature = ComputeGradient (mesh, sol, false, ORDER_2_CENTRAL);

    field->W = BuildWOnBorder (field, mesh, phi, idxsBorder, h0);
    ExtendWToAllDomain (mesh,  idxsBorder, &field->W);

    return field;
}

std::vector<Point*> BuildWOnBorder (Field* field, Mesh* mesh, Vector* phi, std::vector<int>* idxsBorder, double h0)
{
    std::cout << "# Build W vectors on border." << std::endl;

    size_t NumPoints = size_t (mesh->GetNumberOfTotalPoints ());

    std::vector<Point*> W(NumPoints);

    for (size_t i = 0; i < NumPoints; ++i)
        W.at (i) = new Point();

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

            // Gradient de T en p_n
            Point* GradTemperature = field->GradTemperature.at (size_t (p_n->GetGlobalIndex ()));
            Point* Normal = field->Normals.at(size_t (p_n->GetGlobalIndex ()));

            Quantity.at (i) = (*GradTemperature|*Normal);

        }

        // reconnaissance de celui du voisin à l'interieur et du voisin à l'exterieur

        size_t idx_min = 0;
        size_t idx_max = idx_min;

        double phi_min = phi->coeff (listNeigh_p.front ()->GetGlobalIndex ());
        double phi_max = phi_min;

        for (size_t i = 0; i < NumNeigh; ++i)
        {
            Point* p_n = listNeigh_p.at (i);

            double value = phi->coeff (p_n->GetGlobalIndex ());

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

//        double S = -1. / h0 * (Quantity.at (idx_max) - Quantity.at (idx_min));
        double S = -1. / h0 * (Quantity.at (idx_min) - Quantity.at (idx_max));


        Point* Normal_p = field->Normals.at (size_t (idx));

        *W.at (size_t (idx)) = S * *Normal_p;

    }

    std::cout << "\r" << INDENT << "Build W on border is done.       " << std::endl;

    return W;
}

void ExtendWToAllDomain (Mesh* mesh, std::vector<int>* idxsBorder, std::vector<Point*>* W)
{
    std::cout << "# Extend W to all domain." << std::endl;

    // LapW = 0 avec W = W_bordIdxs

    DIM dim = mesh->GetDimension ();
    int NumPoints = mesh->GetNumberOfTotalPoints ();

    Matrix A = Laplacian (mesh);
    //    RemovePeriodicity (mesh, &A);

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

    Vector Wx(NumPoints), Wy(NumPoints), Wz(NumPoints);
    Wx.setZero ();
    Wy.setZero ();
    Wz.setZero ();

    ImposeZeroDirichletExtBorder (mesh, &A, &bx, &by, &bz);

    Wx = Solve (A, bx, IMPLICIT);

    if (dim >= DIM_2D)
        Wy = Solve (A, by, IMPLICIT);

    if (dim == DIM_3D)
        Wz = Solve (A, bz, IMPLICIT);


    for (int i = 0; i < NumPoints; ++i)
        *W->at (size_t (i)) = {Wx.coeff (i), Wy.coeff (i), Wz.coeff (i)};

    std::cout << INDENT << "Extend W is done.       " << std::endl;

    return;
}

std::vector<Point*> ComputeGradient (Mesh* mesh, Vector* vec, bool normalized, ORDERS order)
{
    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hy ();
    double hz = mesh->Get_hz ();

    size_t N = size_t (mesh->GetNumberOfTotalPoints ());

    DIM dim = mesh->GetDimension ();

    std::vector<Point*> R(N);

    for (size_t i = 0; i < N; i++)
        R.at (i) = new Point();

    auto DF = DFStruct();

    std::vector<int> idxDf;
    std::vector<double> coeffDf;

    DFOrderBuild (1, order, &idxDf, &coeffDf);
    size_t SizeDf = idxDf.size ();


    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                Point* p = mesh->GetPoint (i, j, k);
                size_t idx = size_t (p->GetGlobalIndex ());

                Point grad = {0, 0, 0};

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
                        int idx_n = n->GetGlobalIndex ();

                        switch (d)
                        {
                        case DIM_1D:
                            grad.x += coeff * vec->coeff (idx_n) / hx;
                            break;
                        case DIM_2D:
                            grad.y += coeff * vec->coeff (idx_n) / hy;
                            break;
                        case DIM_3D:
                            grad.z += coeff * vec->coeff (idx_n) / hz;
                            break;
                        }
                    }
                }

                double norm = 1.;
                if (normalized)
                    norm = std::sqrt(grad|grad);

                switch (dim)
                {
                case DIM_1D:
                    R.at (idx)->x = grad.x / norm;
                    break;

                case DIM_2D:
                    R.at (idx)->x = grad.x / norm;
                    R.at (idx)->y = grad.y / norm;
                    break;

                case DIM_3D:
                    R.at (idx)->x = grad.x / norm;
                    R.at (idx)->y = grad.y / norm;
                    R.at (idx)->z = grad.z / norm;
                    break;
                }
            }
        }
    }

    Extrapole (mesh, &R);

    return R;
}

Vector ComputeNorm(std::vector<Point*>* vec)
{
    int N = int(vec->size ());
    Vector R(N);
    R.setZero ();

    for (size_t i = 0; i < size_t(N); ++i)
    {
        auto a = *vec->at (i);
        R.coeffRef (int(i)) = std::sqrt(a|a);
    }

    return R;
}

double Compute_dt(Mesh* mesh, Field* field)
{
    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hy ();
    double hz = mesh->Get_hz ();

    DIM dim = mesh->GetDimension ();
    if (dim == DIM_1D)
        hy = hz = 1.;
    if (dim == DIM_2D)
        hz = 1.;

    // dt = min (dt_h, dt_l)
    // dt_h = 0.5 dx or 0.5 dx^2
    // dt_l <= 0.5 / (w1/dx + w2/dy + w3/dz)

    double dt_l = 0.;

    double q = 0.;

    double dt_h = 0.5 * hx;
    //    double dt_h = 0.5 * hx * hx;

    size_t N = field->W.size ();
    for (size_t i = 0; i < N; ++i)
    {
        Point* w = field->W.at (i);
        q = std::max (q, w->x / hx + w->y / hy + w->z / hz);
    }

    dt_l = 0.5 / q;

    return std::min(dt_h, dt_l);
//    return dt_h;


}


