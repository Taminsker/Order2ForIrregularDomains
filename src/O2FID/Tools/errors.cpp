#include "errors.h"

double GetErrorL1 (Mesh * mesh,
                   const Vector &u_ana,
                   const Vector &u_num)
{
    double error_L1;
    double hx = (std::abs(mesh->Get_hx()) < 1e-8) ? 1 : mesh->Get_hx();
    double hy = (std::abs(mesh->Get_hy()) < 1e-8) ? 1 : mesh->Get_hx();
    double hz = (std::abs(mesh->Get_hz()) < 1e-8) ? 1 : mesh->Get_hx();
    Vector u_abs = (u_ana - u_num).cwiseAbs();
    error_L1 = u_abs.sum() * hx * hy * hz;

    return error_L1;
}

double GetErrorLinf (Mesh * mesh,
                     const Vector &u_ana,
                     const Vector &u_num)
{
    Vector u_abs = (u_ana - u_num).cwiseAbs();
    double error_Linf = u_abs.maxCoeff();

    return error_Linf;
}

Vector GetErrorAbs (Mesh * mesh,
                    const Vector &u_ana,
                    const Vector &u_num)
{

}

double GetErrorRela (Mesh * mesh,
                     const Vector &u_ana,
                     const Vector &u_num)
{

}
