#include "solver.h"

#include <Eigen/IterativeLinearSolvers> // pour le gradient conjugué

Vector Solve (const Matrix &A, const Vector &b, TYPE type)
{
    // Verification des tailles de A et b
    if (b.size() != A.cols())
    {
        std::cerr << "Les matries A et b ne sont pas de la même taille !" << std::endl;
    }

    // Déclaration du vecteur solution.
    int size = b.size();
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
            Eigen::ConjugateGradient<Matrix, Eigen::Lower|Eigen::Upper> cg; // déclare la fonction gradient conjugé
            cg.setMaxIterations(int(log(A.size()))); // De base à N.
            cg.setTolerance(1e-8); // De base à la precision de la machine.
            cg.compute(A); // Calcule la matrice A
            sol = cg.solve(b); // Solution du problème
            std::cout << "#iterations:" << cg.iterations() << std::endl;
            std::cout << "erreur estimée:" << cg.error() << std::endl;
            sol = cg.solve(b); // met à jour la solution finale
        }

        if (!A.isApprox(A.adjoint())) // Si non matrice hermitienne on utilise le bi-gradient conjugé
        {
            Eigen::BiCGSTAB<Matrix> solver;
            solver.compute(A);
            solver.setMaxIterations(int(log(A.size()))); // De base à N.
            solver.setTolerance(1e-8); // De base à la precision de la machine.
            sol = solver.solve(b);
            std::cout << "#iterations:     " << solver.iterations() << std::endl;
            std::cout << "estimated error: " << solver.error()      << std::endl;
            /* ... update b ... */
            sol = solver.solve(b); // solve again
        }

    }

    return sol;
}
