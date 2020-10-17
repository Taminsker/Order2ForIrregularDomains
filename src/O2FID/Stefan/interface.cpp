/** @file interface.cpp */

#include "interface.h"


Vector IteratePhi (Mesh* mesh, std::vector<Point*>* W, double dt, Vector* phi_n, TVD tvd)
{

    // Solve phi_t + W \cdot \Grad phi = 0
    Vector phi, u_0, u_1, u_2, u_3;
    u_0 = *phi_n;

    dt = 0.1 * mesh->Get_hx ();

    switch (tvd)
    {
    case TVD_RK2:
        // TVD scheme RK2
        u_1 = u_0 + dt * Weno (mesh, &u_0, W);
        phi = 1. / 2. * u_0 + 1. / 2. * u_1 + dt / 2. *  Weno (mesh, &u_1, W);
        break;

    case TVD_RK3:
        // TVD scheme RK3
        u_1 = u_0 + dt * Weno (mesh, &u_0, W);
        u_2 = 3./4. * u_0 + 1./ 4. * u_1 + 1./ 4. * dt * Weno (mesh, &u_1, W);
        phi = 1./ 3. * u_0 + 2. / 3. * u_2 + 2./ 3. * dt * Weno (mesh, &u_2, W);
        break;

    case TVD_RK4:
        // TVD scheme RK4
        u_1 = u_0 + dt / 2. * Weno (mesh, &u_0, W);
        u_2 = u_1 + dt / 2. * (-1. * Weno(mesh, &u_0, W) + Weno(mesh, &u_1, W));
        u_3 = u_2 + dt / 2. * (-1. *Weno(mesh, &u_1, W) + 2. * Weno(mesh, &u_2, W));
        phi = u_3 + dt / 6. * (Weno(mesh, &u_0, W) + 2. * Weno(mesh, &u_1, W) - 4. * Weno(mesh, &u_2, W) + Weno(mesh, &u_3, W));
        break;
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


                    for (int fl : {-1, 1})
                    {
                        double w_value = 0.;
                        if (fl == -1)
                            w_value = wp;
                        else
                            w_value = wm;

                        for (int inc : {0, 1})
                        {
                            // flux croissant fl == -1, i+1/2 si inc == 0
                            // flux croissant fl == -1, i-1/2 si inc == 1
                            // flux decroissant fl == 1, i+1/2 si inc == 0
                            // flux decroissant fl == 1, i-1/2 si inc == 1

                            int I = i - inc * incI;
                            int J = j - inc * incJ;
                            int K = k - inc * incK;

                            Point *P2, *P1, *ID, *M1, *M2;

                            if (fl == -1)
                            {
                                P2 = mesh->GetPoint (I + 2*incI, J + 2*incJ, K + 2*incK);
                                P1 = mesh->GetPoint (I + incI, J + incJ, K + incK);
                                ID = mesh->GetPoint (I, J, K);
                                M1 = mesh->GetPoint (I - incI, J - incJ, K - incK);
                                M2 = mesh->GetPoint (I - 2*incI, J - 2*incJ, K - 2*incK);
                            } else
                            {
                                P2 = mesh->GetPoint (I - incI, J - incJ, K - incK);
                                P1 = mesh->GetPoint (I, J, k);
                                ID = mesh->GetPoint (I + incI, J + incJ, K + incK);
                                M1 = mesh->GetPoint (I + 2*incI, J + 2*incJ, K + 2*incK);
                                M2 = mesh->GetPoint (I + 3*incI, J + 3*incJ, K + 3*incK);
                            }

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

                            if (inc == 0)
                                *coeff += fl * 1./ h * w_value * flux;
                            else
                                *coeff -= fl * 1./ h * w_value * flux;
                        } // Fin inc
                    } // Fin fl
                } // Fin dim
            } // Fin i
        } // Fin j
    } // Fin k

    std::cout << "\r" << INDENT << "End WENO                " << std::endl;

    std::cout << std::endl;

    return R;
}


Vector ReInitPhi (Mesh* mesh, Vector* phi, int itMax, TVD tvd)
{
    std::cout << "# Reinitialisation of phi" << std::endl;

    Vector u_0, u_1, u_2, u_3;
    u_0 = *phi;

    double dt = 0.5 * mesh->Get_hx ();
    //    int N = mesh->GetNumberOfTotalPoints ();
    //    int G = mesh->GetNumberOfCartesianPoints ();

    for (int i = 0; i < itMax; ++i)
    {
        std::cout << "\r i = " << i << std::flush;
        switch (tvd)
        {
        case TVD_RK2:
            u_1 = u_0 + dt * Hamiltonien (mesh, &u_0, phi);
            u_0 = 1. / 2. * u_0 + 1. / 2. * u_1 + dt / 2. * Hamiltonien (mesh, &u_1, phi);
            break;

        case TVD_RK3:
            u_1 = u_0 + dt * Hamiltonien (mesh, &u_0, phi);
            u_2 = 3./4. * u_0 + 1./ 4. * u_1 + 1./ 4. * dt * Hamiltonien (mesh, &u_1, phi);
            u_0 = 1./ 3. * u_0 + 2. / 3. * u_2 + 2./ 3. * dt * Hamiltonien (mesh, &u_2, phi);
            break;

        case TVD_RK4:
            u_1 = u_0 + dt / 2. * Hamiltonien (mesh, &u_0, phi);
            u_2 = u_1 + dt / 2. * (-1. * Hamiltonien (mesh, &u_0, phi) + Hamiltonien (mesh, &u_1, phi));
            u_3 = u_2 + dt / 2. * (-1. * Hamiltonien (mesh, &u_1, phi) + 2. * Hamiltonien (mesh, &u_2, phi));
            u_0 = u_3 + dt / 6. * (Hamiltonien (mesh, &u_0, phi) + 2. * Hamiltonien (mesh, &u_1, phi) - 4. * Hamiltonien (mesh, &u_2, phi) + Hamiltonien (mesh, &u_3, phi));
            break;
        }
    }

    std::cout << std::endl;

    return u_0;
}

Vector Hamiltonien (Mesh* mesh, Vector* phi, Vector* phi_0)
{
    int N = mesh->GetNumberOfTotalPoints ();

    int Nx = mesh->Get_Nx ();
    int Ny = mesh->Get_Ny ();
    int Nz = mesh->Get_Nz ();

    double hx = mesh->Get_hx ();
    double hy = mesh->Get_hy ();
    double hz = mesh->Get_hz ();

    DIM dim = mesh->GetDimension ();

    double eps = hx;

    Vector SPhi (N);
    for (int idx = 0; idx < N; ++idx)
    {
        double *coeff = &phi_0->coeffRef (idx);
        SPhi.coeffRef (idx) = *coeff / std::sqrt(*coeff * *coeff + eps * eps);
    }

    Vector R(N);
    R.setZero ();

    for (int k = 0; k < Nz; ++k)
    {
        for (int j = 0; j < Ny; ++j)
        {
            for (int i = 0; i < Nx; ++i)
            {
                int idx = mesh->GetPoint (i, j, k)->GetGlobalIndex ();
                double S0 = SPhi.coeff (idx);
                double P0 = phi_0->coeff(idx);
                Point D = {0., 0., 0.};

                for (DIM d : {DIM_1D, DIM_2D, DIM_3D})
                {
                    if (d > dim)
                        break;

                    int I = int(d == DIM_1D);
                    int J = int(d == DIM_2D);
                    int K = int(d == DIM_3D);

                    double D_p = 0.;
                    double D_m = 0.;
                    double h = 0.;
                    double* quantity = nullptr;

                    switch (d)
                    {
                    case DIM_1D:
                        h = hx;
                        quantity = &D.x;
                        break;
                    case DIM_2D:
                        h = hy;
                        quantity = &D.y;
                        break;
                    case DIM_3D:
                        h = hz;
                        quantity = &D.z;
                        break;
                    }

                    if (I+J+K != 1)
                        std::cout << "Error....." << std::endl;

                    double M3 = phi->coeff (mesh->GetPoint (i - 3 * I, j - 3 * J, k - 3 * K)->GetGlobalIndex ());
                    double M2 = phi->coeff (mesh->GetPoint (i - 2 * I, j - 2 * J, k - 2 * K)->GetGlobalIndex ());
                    double M1 = phi->coeff (mesh->GetPoint (i - 1 * I, j - 1 * J, k - 1 * K)->GetGlobalIndex ());
                    double ID = phi->coeff (mesh->GetPoint (i        , j        , k        )->GetGlobalIndex ());
                    double P1 = phi->coeff (mesh->GetPoint (i + 1 * I, j + 1 * J, k + 1 * K)->GetGlobalIndex ());
                    double P2 = phi->coeff (mesh->GetPoint (i + 2 * I, j + 2 * J, k + 2 * K)->GetGlobalIndex ());
                    double P3 = phi->coeff (mesh->GetPoint (i + 3 * I, j + 3 * J, k + 3 * K)->GetGlobalIndex ());

                    double part_1 = 0.;

                    part_1 += -1. * (M1 - M2) / h;
                    part_1 +=  7. * (ID - M1) / h;
                    part_1 +=  7. * (P1 - ID) / h;
                    part_1 += -1. * (P2 - P1) / h;

                    part_1 *= (1./ 12.);

                    // D^+ phi _ijk
                    double part_2 = 0.;

                    double c_a = (P3 - 2. * P2 + P1) / h;
                    double c_b = (P2 - 2. * P1 + ID) / h;
                    double c_c = (P1 - 2. * ID + M1) / h;
                    double c_d = (ID - 2. * M1 + M2) / h;

                    double IS_0 = 13. * (c_a - c_b) * (c_a - c_b) + 3. * (c_a - 3. * c_b) * (c_a - 3. * c_b);
                    double IS_1 = 13. * (c_b - c_c) * (c_b - c_c) + 3. * (c_b + c_c) * (c_b + c_c);
                    double IS_2 = 13. * (c_c - c_d) * (c_c - c_d) + 3. * (3. * c_c - c_d) * (3. * c_c - c_d);

                    double alpha_0 = 1. / ((eps + IS_0) * (eps + IS_0));
                    double alpha_1 = 6. / ((eps + IS_1) * (eps + IS_1));
                    double alpha_2 = 3. / ((eps + IS_2) * (eps + IS_2));

                    double w_0 = alpha_0 / (alpha_0 + alpha_1 + alpha_2);
                    double w_2 = alpha_2 / (alpha_0 + alpha_1 + alpha_2);

                    part_2 = (1. / 3.) * w_0 * (c_a - 2. * c_b + c_c);
                    part_2 += (1. / 6.) * (w_2 - 0.5) * (c_b - 2. * c_c + c_d);

                    D_p = part_1 + part_2;

                    // D^- phi _ijk
                    part_2 = 0.;

                    c_a = (M1 - 2. * M2 + M3) / h;
                    c_b = (ID - 2. * M1 + M2) / h;
                    c_c = (P1 - 2. * ID + M1) / h;
                    c_d = (P2 - 2. * P1 + ID) / h;

                    IS_0 = 13. * (c_a - c_b) * (c_a - c_b) + 3. * (c_a - 3. * c_b) * (c_a - 3. * c_b);
                    IS_1 = 13. * (c_b - c_c) * (c_b - c_c) + 3. * (c_b + c_c) * (c_b + c_c);
                    IS_2 = 13. * (c_c - c_d) * (c_c - c_d) + 3. * (3. * c_c - c_d) * (3. * c_c - c_d);

                    alpha_0 = 1. / ((eps + IS_0) * (eps + IS_0));
                    alpha_1 = 6. / ((eps + IS_1) * (eps + IS_1));
                    alpha_2 = 3. / ((eps + IS_2) * (eps + IS_2));

                    w_0 = alpha_0 / (alpha_0 + alpha_1 + alpha_2);
                    w_2 = alpha_2 / (alpha_0 + alpha_1 + alpha_2);

                    part_2 = (1. / 3.) * w_0 * (c_a - 2. * c_b + c_c);
                    part_2 += (1. / 6.) * (w_2 - 0.5) * (c_b - 2. * c_c + c_d);

                    D_m = part_1 - part_2;

                    if (P0 >= 0)
                        *quantity = std::max(-1. * std::min (D_p, 0.), std::max (D_m, 0.));
                    else
                        *quantity = std::max(std::max (D_p, 0.), -1. * std::min (D_m, 0.));
                }

                R.coeffRef (idx) = S0 * (1. - std::sqrt(D|D));
            }
        }
    }

    return R;
}
