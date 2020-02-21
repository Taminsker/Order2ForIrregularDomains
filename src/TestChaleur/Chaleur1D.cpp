#include "Chaleur1D.h"

Eigen::VectorXd mesh (double a, double b, int n)
{
    Eigen::VectorXd mesh (n);
    double h = (b-a)/double (n + 1);

    for (int i = 0; i < n; i++)
    {
        mesh [i] = a + (i + 1) * h;
    }

    return mesh;
}

Eigen::SparseMatrix<double> matrix (double a, double b, int n, double dt,
                                    double D)
{
    double dx = (b-a)/double (n - 1);
    double alpha = -(D * dt) / (dx * dx);
    double beta = 1 - 2 * alpha;

    /* Remplissage matrice. */
    Eigen::SparseMatrix<double> matrix (n, n);
    matrix.reserve (3 * n);

    for (int i = 0; i < n; i++)
    {
        matrix.insert (i, i) = beta;

        if ((i - 1) >= 0)
        {
            matrix.insert (i, i - 1) = alpha;
        }
        if ((i + 1) < n)
        {
            matrix.insert (i, i + 1) = alpha;
        }
    }

    return matrix;
}

/* Création de la source
 * entrée : maillage de l'espace
  * sortie : vecteur contenant les sources sur tout les points du maillage. */
Eigen::VectorXd source (Eigen::VectorXd mesh)
{
    return 0 * mesh; // Pour le moment une source nulle
}

void solve (Eigen::SparseMatrix<double> matrix, Eigen::VectorXd mesh, double dt
            )
{

}
