#include "interface.h"


Vector IteratePhi (Mesh* mesh, std::vector<Point*>* W, double dt, Vector* phi_n, Vector* phi_n_1, Vector* phi_n_2)
{

    // Solve phi_t + W \cdot \Grad phi = 0

    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hy ();
    double hz = mesh->Get_hz ();

    int N = mesh->GetNumberOfCartesianPoints ();
    DIM dim = mesh->GetDimension ();

    // Partie Spatiale
    Matrix As(N, N);

    auto DF = DFStruct();

    auto idxDf = DF.Derivative_1.Central.Order8.idxs;
    auto coeffDf = DF.Derivative_1.Central.Order8.coeffs;

    if (true)
    {
        size_t SizeDf = idxDf.size ();

        As.reserve (Eigen::VectorXi::Constant(N, int(3 * idxDf.size ())));

        for (int k = 0; k < Nz; ++k)
        {
            for (int j = 0; j < Ny; ++j)
            {
                for (int i = 0; i < Nx; ++i)
                {
                    int idxGlobal = mesh->GetGlobalIndexOfPoint (i, j, k);
                    Point* w = W->at (size_t (idxGlobal));

                    double x = w->x * 1. / hx;
                    double y = w->y * 1. / hy;
                    double z = w->z * 1. / hz;

                    for (size_t s = 0; s < SizeDf; ++s)
                    {
                        double coeff = coeffDf.at (s);
                        int emp = idxDf.at (s);

                        int idx = mesh->GetGlobalIndexOfPoint (i + emp, j, k);
                        As.coeffRef (idxGlobal, idx) = x * coeff;

                        if (dim == DIM_2D || dim == DIM_3D)
                        {
                            idx = mesh->GetGlobalIndexOfPoint (i, j + emp, k);
                            As.coeffRef (idxGlobal, idx) = y * coeff;

                            if (dim == DIM_3D)
                            {
                                idx = mesh->GetGlobalIndexOfPoint (i, j, k + emp);
                                As.coeffRef (idxGlobal, idx) = z * coeff;
                            }
                        }
                    }
                }
            }
        }
    }

    Vector phi;

    if (false)
    {
        Vector b (N);
        b.setZero ();

        // Partie temporelle
        Matrix At(N, N);
        At.setIdentity ();

        // Ordre 1
        //    At *= 1. / m_dt;
        //    b += At * *m_chain->front ();

        //Ordre 2
        //    At *= (3./ 2.)  / m_dt;
        //    b -= - 2.       / m_dt * *m_chain->back ();
        //    b -=  (1./ 2.)  / m_dt * *m_chain->at (m_chain->size ()-2);

        //Ordre 3
        At *= (11./ 6.) / dt;
        b -= - 3.       / dt * *phi_n;
        b -= (3./ 2.)   / dt * *phi_n_1;
        b -= (-1./3.)   / dt * *phi_n_2;

        //    Matrix A = At + 0.5 * As;
        //        b -= (0.5 * As) * *m_chain->back ();
        //    b += At * *m_chain->back ();

        Matrix A = At + As;

        phi = Solve (A, b, IMPLICIT);
    }

    if (true)
    {
        // TVD scheme
        Vector u_n = *phi_n;

        Vector u_1 = u_n - dt * As * u_n;
        Vector u_2 = 3./4. * u_n + 1./ 4. * u_1 - 1./ 4. * dt * As * u_1;
        phi = 1./ 3. * u_n + 2. / 3. * u_2 - 2./ 3. * dt * As * u_2;
    }

    if (false)
    {
        // TVD scheme
        Vector u_n = *phi_n;

        Vector u_1 = u_n - dt * Weno (mesh, &u_n, W);
        Vector u_2 = 3./4. * u_n + 1./ 4. * u_1 - 1./ 4. * dt * Weno (mesh, &u_1, W);
        phi = 1./ 3. * u_n + 2. / 3. * u_2 - 2./ 3. * dt * Weno (mesh, &u_2, W);
    }

    Extrapole (mesh, &phi);


    std::cout << std::endl;

    return phi;
}

Vector Weno(Mesh* mesh, Vector* u, std::vector<Point*>* W)
{
    std::cout << "# WENO scheme" << std::endl;

    int N = mesh->GetNumberOfCartesianPoints ();

    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hy ();
    double hz = mesh->Get_hz ();

    DIM dim = mesh->GetDimension ();

    Vector R(N);
    R.setZero ();

    double eps = 1e-6;
    double gamma_1 = 1. / 10.;
    double gamma_2 = 3. / 5.;
    double gamma_3 = 3. / 10.;

    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                int idx = mesh->GetGlobalIndexOfPoint (i, j, k);

                //                std::cout << "\r" << INDENT << "Point : (" << i << ", " << j << ", " << k << ") / (" << Nx << ", " << Ny << ", " << Nz << ")            " << std::flush;

                std::cout << "\r" << INDENT << "Point : " << idx << "/" << N << "            " << std::flush;

                Point* w = W->at (size_t (idx));

                double* coeff = &R.coeffRef (idx);

                *coeff = 0.;

                for (DIM d : {DIM_1D, DIM_2D, DIM_3D})
                {
                    if (d > dim)
                        break;

                    double tempCoeff = 0.;

                    for (int fl : {0, -1})
                    {
                        int incI = int(d == DIM_1D);
                        int incJ = int(d == DIM_2D);
                        int incK = int(d == DIM_3D);

                        int I = i + incI * fl;
                        int J = j + incJ * fl;
                        int K = k + incK * fl;

                        // Calcul sur la dimension d

                        int P2 = mesh->GetGlobalIndexOfPoint (I + 2*incI, J + 2*incJ, K + 2*incK);
                        int P1 = mesh->GetGlobalIndexOfPoint (I + incI, J + incJ, K + incK);
                        int ID = mesh->GetGlobalIndexOfPoint (I, J, K);
                        int M1 = mesh->GetGlobalIndexOfPoint (I - incI, J - incJ, K - incK);
                        int M2 = mesh->GetGlobalIndexOfPoint (I - 2*incI, J - 2*incJ, K - 2*incK);

                        double* CP2 = &u->coeffRef (P2);
                        double* CP1 = &u->coeffRef (P1);
                        double* CID = &u->coeffRef (ID);
                        double* CM1 = &u->coeffRef (M1);
                        double* CM2 = &u->coeffRef (M2);


                        double temp = (*CM2 - 2. * *CM1 + *CID);
                        double beta_1 = 13. / 12. * temp * temp;

                        temp = (*CM2 - 4. * *CM1 + 3. * *CID);
                        beta_1 += 1. / 4. * temp * temp;

                        temp = (*CM1 - 2. * *CID + *CP1);
                        double beta_2 = 13. / 12. * temp * temp;

                        temp = (*CM1 - *CP1);
                        beta_2 += 1. / 4. * temp * temp;

                        temp = (*CID - 2. * *CP1 + *CP2);
                        double beta_3 = 13. / 12. * temp * temp;

                        temp = (3. * *CID - 4. * *CP1 + *CP2);
                        beta_3 += 1. / 4. * temp * temp;

                        double omegaB_1 = gamma_1 / ((eps + beta_1) * (eps + beta_1));
                        double omegaB_2 = gamma_2 / ((eps + beta_2) * (eps + beta_2));
                        double omegaB_3 = gamma_3 / ((eps + beta_3) * (eps + beta_3));

                        double sum = omegaB_1 + omegaB_2 + omegaB_3;
                        double omega_1 = omegaB_1 / sum;
                        double omega_2 = omegaB_2 / sum;
                        double omega_3 = omegaB_3 / sum;

                        double f1 = 1. / 3. * *CM2;
                        f1 -= 7./ 6. * *CM1;
                        f1 += 11. / 6. * *CID;

                        double f2 = -1. / 6. * *CM1;
                        f2 += 5./ 6. * *CID;
                        f2 += 1. / 3. * *CP1;

                        double f3 = 1. / 3. * *CID;
                        f3 += 5./ 6. * *CP1;
                        f3 -= 1. / 6. * *CP2;

                        if (fl == 0)
                            tempCoeff += omega_1 * f1 + omega_2 * f2 + omega_3 * f3;
                        else
                            tempCoeff -= omega_1 * f1 + omega_2 * f2 + omega_3 * f3;
                    }

                    switch (d)
                    {
                    case DIM_1D:
                        *coeff += tempCoeff * w->x / hx;
                        break;
                    case DIM_2D:
                        *coeff += tempCoeff * w->y / hy;
                        break;
                    case DIM_3D:
                        *coeff += tempCoeff * w->z / hz;
                        break;
                    }
                }
            }
        }
    }

    std::cout << "\r" << INDENT << "End WENO                " << std::endl;

    std::cout << std::endl;

    return R;
}

void ReInitPhi (Mesh* mesh, Vector* phi, std::vector<int>* idxsBorder, double dt)
{
    std::cout << "# Reinitialisation of phi" << std::endl;

    *idxsBorder = mesh->GetListOfIndexPoints ();
    int N = mesh->GetNumberOfTotalPoints ();
    int G = mesh->GetNumberOfCartesianPoints ();


    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hy ();
    double hz = mesh->Get_hz ();

    DIM dim = mesh->GetDimension ();


    //    std::cout << Phi->transpose () << std::endl;


    double eps = hx;
    if (dim >= DIM_2D)
    {
        eps = std::min(eps, hy);
        if (dim == DIM_3D)
            eps = std::min(eps, hz);
    }

    Vector SPhi (N);
    for (int idx = 0; idx < N; ++idx)
    {
        double *coeff = &phi->coeffRef (idx);
        SPhi.coeffRef (idx) = *coeff / (*coeff * *coeff + eps * eps);
    }

    // Scheme Euler Explicit Order 1 time, Order 2 in space

    Vector Rn;
    Vector Rn_1;


    Extrapole (mesh, phi);

    Rn = *phi;
    Rn_1 = Rn;


    int Limit = 6;
    for (int timer = 1; timer <= Limit; timer++)
    {
        for (int k = 0; k < Nz; ++k)
        {
            for (int j = 0; j < Ny; ++j)
            {
                for (int i = 0; i < Nx; ++i)
                {
                    Point* p = mesh->GetPoint (i, j ,k);

                    int idx = p->GetGlobalIndex ();

                    std::cout << "\r" << INDENT << "(" << timer << "/" << Limit << ") Point : " << idx+1 << "/" << G << "            " << std::flush;

                    Point Grad;
                    Point dist;

                    auto neigh = p->GetListNeighbours ();

                    for (Point* n : neigh)
                    {
                        Point diff = *p - *n;

                        if (diff == Point(diff.x, 0, 0))
                        {
                            if (diff.x < 0 && fabs(diff.x) < 4. * hx)
                                continue;

                            Grad.x += Rn_1.coeff (n->GetGlobalIndex ());
                            dist.x += fabs(diff.x);

                        } else if (diff == Point(0, diff.y, 0) && dim >= DIM_2D)
                        {
                            if (diff.y < 0 && fabs(diff.y) < 4. * hy)
                                continue;

                            Grad.y += Rn_1.coeff (n->GetGlobalIndex ());
                            dist.y += fabs(diff.y);

                        } else if (diff == Point(0, 0, diff.z) && dim == DIM_3D)
                        {
                            if (diff.z < 0 && fabs(diff.z) < 4. * hz)
                                continue;

                            Grad.z += Rn_1.coeff (n->GetGlobalIndex ());
                            dist.z += fabs(diff.z);
                        }
                        else
                        {
                            std::cerr << "ERROR" << std::endl;
                            exit(0);
                        }
                    }

                    double* coeff = &Rn_1.coeffRef (idx);

                    Grad.x = (Grad.x - *coeff) / dist.x;
                    if (dim >= DIM_2D)
                        Grad.y = (Grad.y - *coeff) / dist.y;
                    if (dim == DIM_3D)
                        Grad.z = (Grad.z - *coeff) / dist.z;


                    Rn.coeffRef (idx) = *coeff - dt * SPhi.coeff (idx) * (std::sqrt(Grad|Grad) - 1.);
//                    std::cout << "\nGrad : " << Grad << std::endl;
//                    std::cout << "dist : " << dist << std::endl;

                }
            }
        }

        for (int idx : *idxsBorder)
        {
            Rn.coeffRef (idx) = phi->coeff (idx);
//            std::cout << Phi->coeff (idx) << " " << std::flush;

        }

//        std::cout << std::endl;

        Rn_1 = Rn;
    }

    std::cout << "\r" << INDENT << "End Reinitialisation phi                          " << std::endl;
    std::cout << std::endl;

    *phi = Rn;

//        std::cout << Phi->transpose ()<< std::endl;



    return;
}
