#include "solver.h"

#include <Eigen/IterativeLinearSolvers> // pour le gradient conjugué

Vector Solve (const Matrix &A, const Vector &b, TYPE type)
{
    if (type == IMPLICIT)
    {
        std::cout << "# Solver implicit" << std::endl;
    } else
    {
        std::cout << "# Solver explicit" << std::endl;
    }

    // Verification des tailles de A et b
    if (b.size() != A.cols() && b.transpose ().size() != A.cols())
    {
        std::cerr << INDENT << "ERROR::Solver matrix and vector are not the same size." << std::endl;
        exit(0);
    }

    // Déclaration du vecteur solution.
    int size = int (b.size());
    Vector sol(size);


    // Deux cas de figures : Explicite et Implicite
    if (type == EXPLICIT) // U = Ab : produit matrice vecteur avec U = sol
    {
        sol = A * b;
    }
    if (type == IMPLICIT) // AU = b : gradient conjugé
    {
        if (A.isApprox(A.adjoint())) // Si matrice hermitienne on utilise le gradient conjugé
        {
            std::cout << INDENT << "Symmetric matrix detected." <<  std::endl;

            Eigen::ConjugateGradient<Matrix, Eigen::Lower|Eigen::Upper> cg; // déclare la fonction gradient conjugé
//            cg.setMaxIterations(int(log(A.size()))); // De base à N.
            cg.setTolerance(1e-18); // De base à la precision de la machine.
            cg.compute(A); // Calcule la matrice A
            sol = cg.solve(b); // Solution du problème
            std::cout << INDENT << "Iterations        :" << cg.iterations() << std::endl;
            std::cout << INDENT << "Estimated error   :" << cg.error() << std::endl;
//            sol = cg.solve(b); // met à jour la solution finale
        }
        else // Si non matrice hermitienne on utilise le bi-gradient conjugé
        {
            std::cout << INDENT << "Asymmetric matrix detected." <<  std::endl;

            Eigen::BiCGSTAB<Matrix> solver;
            solver.compute(A);
//            solver.setMaxIterations(int(log(A.size()))); // De base à N.
            solver.setTolerance(1e-18); // De base à la precision de la machine.
            sol = solver.solve(b);
            std::cout << INDENT << "Iterations :        " << solver.iterations() << std::endl;
            std::cout << INDENT << "Estimated error :   " << solver.error()      << std::endl;
            /* ... update b ... */
//            sol = solver.solve(b); // solve again
        }
    }

    std::cout << std::endl;
    return sol;
}
