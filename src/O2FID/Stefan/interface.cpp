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

    auto idxDf = DF.Derivative_1.Central.Order2.idxs;
    auto coeffDf = DF.Derivative_1.Central.Order2.coeffs;

    if (false)
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

        RemovePeriodicity (mesh, &As);
    }

    Vector phi;

    if (false)
    {
        Vector b (N);
        b.setZero ();

        // Partie temporelle
        Matrix At(N, N);
        At.setIdentity ();

        //         Ordre 1
        At *= 1. / dt;
        b += At * *phi_n;

        //Ordre 2
        //    At *= (3./ 2.)  / m_dt;
        //    b -= - 2.       / m_dt * *m_chain->back ();
        //    b -=  (1./ 2.)  / m_dt * *m_chain->at (m_chain->size ()-2);

        //        //Ordre 3
        //        At *= (11./ 6.) / dt;
        //        b -= - 3.       / dt * *phi_n;
        //        b -= (3./ 2.)   / dt * *phi_n_1;
        //        b -= (-1./3.)   / dt * *phi_n_2;

        //    Matrix A = At + 0.5 * As;
        //        b -= (0.5 * As) * *m_chain->back ();
        //    b += At * *m_chain->back ();

        Matrix A = At + As;

        phi = Solve (A, b, IMPLICIT);
    }

    if (false)
    {
        // TVD scheme
        Vector u_n = *phi_n;

        Vector u_1 = u_n - dt * As * u_n;
        Vector u_2 = 3./4. * u_n + 1./ 4. * u_1 - 1./ 4. * dt * As * u_1;
        phi = 1./ 3. * u_n + 2. / 3. * u_2 - 2./ 3. * dt * As * u_2;
    }

    if (true)
    {
        // TVD scheme
        Vector u_n = *phi_n;

        Vector u_1 = u_n + dt * Weno (mesh, &u_n, W);
        Vector u_2 = 3./4. * u_n + 1./ 4. * u_1 + 1./ 4. * dt * Weno (mesh, &u_1, W);
        phi = 1./ 3. * u_n + 2. / 3. * u_2 + 2./ 3. * dt * Weno (mesh, &u_2, W);
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

    double maxima = std::max(hx, std::max(hy, hz)) + 1e-10;

    DIM dim = mesh->GetDimension ();

    Vector R(N);
    R.setZero ();

    double eps = 1e-6;
    double gamma_1 = 1. / 10.;
    double gamma_2 = 3. / 5.;
    double gamma_3 = 3. / 10.;

//    double ZERO = 0.;

    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                int idx = mesh->GetGlobalIndexOfPoint (i, j, k);

//                Point* p = mesh->GetPoint (idx);

                //                std::cout << "\r" << INDENT << "Point : (" << i << ", " << j << ", " << k << ") / (" << Nx << ", " << Ny << ", " << Nz << ")            " << std::flush;

                std::cout << "\r" << INDENT << "Point : " << idx << "/" << N << "            " << std::flush;

                Point* w = W->at (size_t (idx));

                double* coeff = &R.coeffRef (idx);

                *coeff = 0.;

                for (DIM d : {DIM_1D, DIM_2D, DIM_3D})
                {
                    if (d > dim)
                        break;

                    double wp = 0.;
                    double wm = 0.;
                    double h = 0.;

                    switch (d)
                    {
                    case DIM_1D:
                        h = hx;
                        wp = (fabs(w->x) + w->x) / 2.;
                        wm = (fabs(w->x) - w->x) / 2.;
                        break;
                    case DIM_2D:
                        wp = (fabs(w->y) + w->y) / 2.;
                        wm = (fabs(w->y) - w->y) / 2.;
                        h = hy;
                        break;
                    case DIM_3D:
                        wp = (fabs(w->z) + w->z) / 2.;
                        wm = (fabs(w->z) - w->z) / 2.;
                        h = hz;
                        break;
                    }

                    int incI = int(d == DIM_1D);
                    int incJ = int(d == DIM_2D);
                    int incK = int(d == DIM_3D);


                    for (int fl : {0, 1})
                    {
                        // Partie flux croissant, i+1/2

                        int I = i - fl * incI;
                        int J = j - fl * incJ;
                        int K = k - fl * incK;

                        Point* P2 = mesh->GetPoint (I + 2*incI, J + 2*incJ, K + 2*incK);
                        Point* P1 = mesh->GetPoint (I + incI, J + incJ, K + incK);
                        Point* ID = mesh->GetPoint (I, J, K);
                        Point* M1 = mesh->GetPoint (I - incI, J - incJ, K - incK);
                        Point* M2 = mesh->GetPoint (I - 2*incI, J - 2*incJ, K - 2*incK);

                        double CP2 = u->coeffRef (P2->GetGlobalIndex ());
                        double CP1 = u->coeffRef (P1->GetGlobalIndex ());
                        double CID = u->coeffRef (ID->GetGlobalIndex ());
                        double CM1 = u->coeffRef (M1->GetGlobalIndex ());
                        double CM2 = u->coeffRef (M2->GetGlobalIndex ());


                        double temp = (CM2 - 2. * CM1 + CID);
                        double beta_1 = 13. / 12. * temp * temp;

                        temp = (CM2 - 4. * CM1 + 3. * CID);
                        beta_1 += 1. / 4. * temp * temp;

                        temp = (CM1 - 2. * CID + CP1);
                        double beta_2 = 13. / 12. * temp * temp;

                        temp = (CM1 - CP1);
                        beta_2 += 1. / 4. * temp * temp;

                        temp = (CID - 2. * CP1 + CP2);
                        double beta_3 = 13. / 12. * temp * temp;

                        temp = (3. * CID - 4. * CP1 + CP2);
                        beta_3 += 1. / 4. * temp * temp;

                        double omegaB_1 = gamma_1 / ((eps + beta_1) * (eps + beta_1));
                        double omegaB_2 = gamma_2 / ((eps + beta_2) * (eps + beta_2));
                        double omegaB_3 = gamma_3 / ((eps + beta_3) * (eps + beta_3));

                        double sum = omegaB_1 + omegaB_2 + omegaB_3;
                        double omega_1 = omegaB_1 / sum;
                        double omega_2 = omegaB_2 / sum;
                        double omega_3 = omegaB_3 / sum;

                        double f1 = 1. / 6. * (2.  * CM2 - 7. * CM1 + 11. * CID);
                        double f2 = 1. / 6. * (-1. * CM1 + 5. * CID + 2.  * CP1);
                        double f3 = 1. / 6. * (2.  * CID + 5. * CP1 - 1.  * CP2);

                        double flux = omega_1 * f1 + omega_2 * f2 + omega_3 * f3;

                        if (fl == 0)
                            *coeff -= 1./ h * wp * flux;
                        else
                            *coeff += 1./ h * wp * flux;

                        // Partie flux decroissants


                        P2 = mesh->GetPoint (I - incI, J - incJ, K - incK);
                        P1 = mesh->GetPoint (I, J, k);
                        ID = mesh->GetPoint (I + incI, J + incJ, K + incK);
                        M1 = mesh->GetPoint (I + 2*incI, J + 2*incJ, K + 2*incK);
                        M2 = mesh->GetPoint (I + 3*incI, J + 3*incJ, K + 3*incK);

                        CP2 = u->coeffRef (P2->GetGlobalIndex ());
                        CP1 = u->coeffRef (P1->GetGlobalIndex ());
                        CID = u->coeffRef (ID->GetGlobalIndex ());
                        CM1 = u->coeffRef (M1->GetGlobalIndex ());
                        CM2 = u->coeffRef (M2->GetGlobalIndex ());


                        temp = (CM2 - 2. * CM1 + CID);
                        beta_1 = 13. / 12. * temp * temp;

                        temp = (CM2 - 4. * CM1 + 3. * CID);
                        beta_1 += 1. / 4. * temp * temp;

                        temp = (CM1 - 2. * CID + CP1);
                        beta_2 = 13. / 12. * temp * temp;

                        temp = (CM1 - CP1);
                        beta_2 += 1. / 4. * temp * temp;

                        temp = (CID - 2. * CP1 + CP2);
                        beta_3 = 13. / 12. * temp * temp;

                        temp = (3. * CID - 4. * CP1 + CP2);
                        beta_3 += 1. / 4. * temp * temp;

                        omegaB_1 = gamma_1 / ((eps + beta_1) * (eps + beta_1));
                        omegaB_2 = gamma_2 / ((eps + beta_2) * (eps + beta_2));
                        omegaB_3 = gamma_3 / ((eps + beta_3) * (eps + beta_3));

                        sum = omegaB_1 + omegaB_2 + omegaB_3;
                        omega_1 = omegaB_1 / sum;
                        omega_2 = omegaB_2 / sum;
                        omega_3 = omegaB_3 / sum;

                        f1 = 1. / 6. * (2.  * CM2 - 7. * CM1 + 11. * CID);
                        f2 = 1. / 6. * (-1. * CM1 + 5. * CID + 2.  * CP1);
                        f3 = 1. / 6. * (2.  * CID + 5. * CP1 - 1.  * CP2);

                        flux = omega_1 * f1 + omega_2 * f2 + omega_3 * f3;

                        if (fl == 0)
                            *coeff += 1./ h * wm * flux;
                        else
                            *coeff -= 1./ h * wm * flux;
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


    int Limit = 3;
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

Vector GradNorm(Mesh* mesh, Vector* u)
{

    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hy ();
    double hz = mesh->Get_hz ();

    int N = mesh->GetNumberOfTotalPoints ();
    DIM dim = mesh->GetDimension ();

    Vector R(N);
    R.setZero ();

    auto DF = DFStruct();

    auto idxDf = DF.Derivative_1.Central.Order8.idxs;
    auto coeffDf = DF.Derivative_1.Central.Order8.coeffs;

    size_t SizeDf = idxDf.size ();

    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                int idxGlobal = mesh->GetGlobalIndexOfPoint (i, j, k);

                Point Grad;

                for (size_t s = 0; s < SizeDf; ++s)
                {
                    double coeff = coeffDf.at (s);
                    int emp = idxDf.at (s);

                    int idx = mesh->GetGlobalIndexOfPoint (i + emp, j, k);
                    double* value = &u->coeffRef (idx);

                    Grad.x += coeff * *value / hx;

                    if (dim == DIM_2D || dim == DIM_3D)
                    {
                        idx = mesh->GetGlobalIndexOfPoint (i, j + emp, k);
                        Grad.y += coeff * *value / hy;

                        if (dim == DIM_3D)
                        {
                            Grad.z += coeff * *value / hz;
                        }
                    }
                }

                R.coeffRef (idxGlobal) = std::sqrt(Grad|Grad);
            }
        }
    }

    Extrapole (mesh, &R);

    return R;
}


//// Partie flux croissant, i+1/2
//{
//    int I = i;
//    int J = j;
//    int K = k;

//    Point* P2 = mesh->GetPoint (I + 2*incI, J + 2*incJ, K + 2*incK);
//    Point* P1 = mesh->GetPoint (I + incI, J + incJ, K + incK);
//    Point* ID = mesh->GetPoint (I, J, K);
//    Point* M1 = mesh->GetPoint (I - incI, J - incJ, K - incK);
//    Point* M2 = mesh->GetPoint (I - 2*incI, J - 2*incJ, K - 2*incK);

//    double CP2 = u->coeffRef (P2->GetGlobalIndex ());
//    double CP1 = u->coeffRef (P1->GetGlobalIndex ());
//    double CID = u->coeffRef (ID->GetGlobalIndex ());
//    double CM1 = u->coeffRef (M1->GetGlobalIndex ());
//    double CM2 = u->coeffRef (M2->GetGlobalIndex ());


//    double temp = (CM2 - 2. * CM1 + CID);
//    double beta_1 = 13. / 12. * temp * temp;

//    temp = (CM2 - 4. * CM1 + 3. * CID);
//    beta_1 += 1. / 4. * temp * temp;

//    temp = (CM1 - 2. * CID + CP1);
//    double beta_2 = 13. / 12. * temp * temp;

//    temp = (CM1 - CP1);
//    beta_2 += 1. / 4. * temp * temp;

//    temp = (CID - 2. * CP1 + CP2);
//    double beta_3 = 13. / 12. * temp * temp;

//    temp = (3. * CID - 4. * CP1 + CP2);
//    beta_3 += 1. / 4. * temp * temp;

//    double omegaB_1 = gamma_1 / ((eps + beta_1) * (eps + beta_1));
//    double omegaB_2 = gamma_2 / ((eps + beta_2) * (eps + beta_2));
//    double omegaB_3 = gamma_3 / ((eps + beta_3) * (eps + beta_3));

//    double sum = omegaB_1 + omegaB_2 + omegaB_3;
//    double omega_1 = omegaB_1 / sum;
//    double omega_2 = omegaB_2 / sum;
//    double omega_3 = omegaB_3 / sum;

//    double f1 = 1. / 6. * (2.  * CM2 - 7. * CM1 + 11. * CID);
//    double f2 = 1. / 6. * (-1. * CM1 + 5. * CID + 2.  * CP1);
//    double f3 = 1. / 6. * (2.  * CID + 5. * CP1 - 1.  * CP2);

//    double flux = omega_1 * f1 + omega_2 * f2 + omega_3 * f3;

//    *coeff -= 1./ h * wp * flux;
//}

//// Partie flux croissant, i-1/2
//{
//    int I = i - incI;
//    int J = j - incJ;
//    int K = k - incK;

//    Point* P2 = mesh->GetPoint (I + 2*incI, J + 2*incJ, K + 2*incK);
//    Point* P1 = mesh->GetPoint (I + incI, J + incJ, K + incK);
//    Point* ID = mesh->GetPoint (I, J, K);
//    Point* M1 = mesh->GetPoint (I - incI, J - incJ, K - incK);
//    Point* M2 = mesh->GetPoint (I - 2*incI, J - 2*incJ, K - 2*incK);

//    double CP2 = u->coeffRef (P2->GetGlobalIndex ());
//    double CP1 = u->coeffRef (P1->GetGlobalIndex ());
//    double CID = u->coeffRef (ID->GetGlobalIndex ());
//    double CM1 = u->coeffRef (M1->GetGlobalIndex ());
//    double CM2 = u->coeffRef (M2->GetGlobalIndex ());


//    double temp = (CM2 - 2. * CM1 + CID);
//    double beta_1 = 13. / 12. * temp * temp;

//    temp = (CM2 - 4. * CM1 + 3. * CID);
//    beta_1 += 1. / 4. * temp * temp;

//    temp = (CM1 - 2. * CID + CP1);
//    double beta_2 = 13. / 12. * temp * temp;

//    temp = (CM1 - CP1);
//    beta_2 += 1. / 4. * temp * temp;

//    temp = (CID - 2. * CP1 + CP2);
//    double beta_3 = 13. / 12. * temp * temp;

//    temp = (3. * CID - 4. * CP1 + CP2);
//    beta_3 += 1. / 4. * temp * temp;

//    double omegaB_1 = gamma_1 / ((eps + beta_1) * (eps + beta_1));
//    double omegaB_2 = gamma_2 / ((eps + beta_2) * (eps + beta_2));
//    double omegaB_3 = gamma_3 / ((eps + beta_3) * (eps + beta_3));

//    double sum = omegaB_1 + omegaB_2 + omegaB_3;
//    double omega_1 = omegaB_1 / sum;
//    double omega_2 = omegaB_2 / sum;
//    double omega_3 = omegaB_3 / sum;

//    double f1 = 1. / 6. * (2.  * CM2 - 7. * CM1 + 11. * CID);
//    double f2 = 1. / 6. * (-1. * CM1 + 5. * CID + 2.  * CP1);
//    double f3 = 1. / 6. * (2.  * CID + 5. * CP1 - 1.  * CP2);

//    double flux = omega_1 * f1 + omega_2 * f2 + omega_3 * f3;

//    *coeff += 1./ h * wp * flux;
//}


//// Partie flux decroissant, i+1/2
//{
//    int I = i;
//    int J = j;
//    int K = k;

//    Point* P2 = mesh->GetPoint (I - incI, J - incJ, K - incK);
//    Point* P1 = mesh->GetPoint (I, J, k);
//    Point* ID = mesh->GetPoint (I + incI, J + incJ, K + incK);
//    Point* M1 = mesh->GetPoint (I + 2*incI, J + 2*incJ, K + 2*incK);
//    Point* M2 = mesh->GetPoint (I + 3*incI, J + 3*incJ, K + 3*incK);

//    double CP2 = u->coeffRef (P2->GetGlobalIndex ());
//    double CP1 = u->coeffRef (P1->GetGlobalIndex ());
//    double CID = u->coeffRef (ID->GetGlobalIndex ());
//    double CM1 = u->coeffRef (M1->GetGlobalIndex ());
//    double CM2 = u->coeffRef (M2->GetGlobalIndex ());


//    double temp = (CM2 - 2. * CM1 + CID);
//    double beta_1 = 13. / 12. * temp * temp;

//    temp = (CM2 - 4. * CM1 + 3. * CID);
//    beta_1 += 1. / 4. * temp * temp;

//    temp = (CM1 - 2. * CID + CP1);
//    double beta_2 = 13. / 12. * temp * temp;

//    temp = (CM1 - CP1);
//    beta_2 += 1. / 4. * temp * temp;

//    temp = (CID - 2. * CP1 + CP2);
//    double beta_3 = 13. / 12. * temp * temp;

//    temp = (3. * CID - 4. * CP1 + CP2);
//    beta_3 += 1. / 4. * temp * temp;

//    double omegaB_1 = gamma_1 / ((eps + beta_1) * (eps + beta_1));
//    double omegaB_2 = gamma_2 / ((eps + beta_2) * (eps + beta_2));
//    double omegaB_3 = gamma_3 / ((eps + beta_3) * (eps + beta_3));

//    double sum = omegaB_1 + omegaB_2 + omegaB_3;
//    double omega_1 = omegaB_1 / sum;
//    double omega_2 = omegaB_2 / sum;
//    double omega_3 = omegaB_3 / sum;

//    double f1 = 1. / 6. * (2.  * CM2 - 7. * CM1 + 11. * CID);
//    double f2 = 1. / 6. * (-1. * CM1 + 5. * CID + 2.  * CP1);
//    double f3 = 1. / 6. * (2.  * CID + 5. * CP1 - 1.  * CP2);

//    double flux = omega_1 * f1 + omega_2 * f2 + omega_3 * f3;

//    *coeff += 1./ h * wm * flux;
//}

//// Partie flux decroissant, i-1/2
//{
//    int I = i - incI;
//    int J = j - incJ;
//    int K = k - incK;

//    Point* P2 = mesh->GetPoint (I - incI, J - incJ, K - incK);
//    Point* P1 = mesh->GetPoint (I, J, k);
//    Point* ID = mesh->GetPoint (I + incI, J + incJ, K + incK);
//    Point* M1 = mesh->GetPoint (I + 2*incI, J + 2*incJ, K + 2*incK);
//    Point* M2 = mesh->GetPoint (I + 3*incI, J + 3*incJ, K + 3*incK);

//    double CP2 = u->coeffRef (P2->GetGlobalIndex ());
//    double CP1 = u->coeffRef (P1->GetGlobalIndex ());
//    double CID = u->coeffRef (ID->GetGlobalIndex ());
//    double CM1 = u->coeffRef (M1->GetGlobalIndex ());
//    double CM2 = u->coeffRef (M2->GetGlobalIndex ());


//    double temp = (CM2 - 2. * CM1 + CID);
//    double beta_1 = 13. / 12. * temp * temp;

//    temp = (CM2 - 4. * CM1 + 3. * CID);
//    beta_1 += 1. / 4. * temp * temp;

//    temp = (CM1 - 2. * CID + CP1);
//    double beta_2 = 13. / 12. * temp * temp;

//    temp = (CM1 - CP1);
//    beta_2 += 1. / 4. * temp * temp;

//    temp = (CID - 2. * CP1 + CP2);
//    double beta_3 = 13. / 12. * temp * temp;

//    temp = (3. * CID - 4. * CP1 + CP2);
//    beta_3 += 1. / 4. * temp * temp;

//    double omegaB_1 = gamma_1 / ((eps + beta_1) * (eps + beta_1));
//    double omegaB_2 = gamma_2 / ((eps + beta_2) * (eps + beta_2));
//    double omegaB_3 = gamma_3 / ((eps + beta_3) * (eps + beta_3));

//    double sum = omegaB_1 + omegaB_2 + omegaB_3;
//    double omega_1 = omegaB_1 / sum;
//    double omega_2 = omegaB_2 / sum;
//    double omega_3 = omegaB_3 / sum;

//    double f1 = 1. / 6. * (2.  * CM2 - 7. * CM1 + 11. * CID);
//    double f2 = 1. / 6. * (-1. * CM1 + 5. * CID + 2.  * CP1);
//    double f3 = 1. / 6. * (2.  * CID + 5. * CP1 - 1.  * CP2);

//    double flux = omega_1 * f1 + omega_2 * f2 + omega_3 * f3;

//    *coeff -= 1./ h * wm * flux;
//}
