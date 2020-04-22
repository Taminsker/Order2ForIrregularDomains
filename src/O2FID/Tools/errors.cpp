#include "errors.h"

double GetErrorl1 (Mesh * mesh,
                   const Vector &u_ana,
                   const Vector &u_num)
{
    double error_l1;
    double hx = (std::abs(mesh->Get_hx()) < 1e-8) ? 1 : mesh->Get_hx();
    double hy = (std::abs(mesh->Get_hy()) < 1e-8) ? 1 : mesh->Get_hy();
    double hz = (std::abs(mesh->Get_hz()) < 1e-8) ? 1 : mesh->Get_hz();
    Vector u_abs = (u_ana - u_num).cwiseAbs();
    error_l1 = u_abs.sum() * hx * hy * hz;

    std::cout << "# Error l^1 has been calculated : " << error_l1 << std::endl;

    return error_l1;
}

double GetErrorlinf (Mesh * mesh,
                     const Vector &u_ana,
                     const Vector &u_num)
{
    (void)mesh;

    Vector u_abs = (u_ana - u_num).cwiseAbs();
    double error_linf = u_abs.maxCoeff();

    std::cout << "# Error l^inf has been calculated : " << error_linf << std::endl;

    return error_linf;
}

Vector GetErrorAbs (Mesh * mesh,
                    const Vector &u_ana,
                    const Vector &u_num)
{
    (void)mesh;

    std::cout << "# Error abs has been calculated." << std::endl;
    return (u_ana - u_num).cwiseAbs ();
}

double GetErrorRela (Mesh * mesh,
                     const Vector &u_ana,
                     const Vector &u_num)
{
    (void)mesh;

    double error_rela = ((u_ana - u_num).squaredNorm ()) / u_ana.squaredNorm ();
    std::cout << "# Error rela has been calculated : " << error_rela << std::endl;
    return error_rela;
}
