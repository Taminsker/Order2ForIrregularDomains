#include "solver.h"

#include <Eigen/IterativeLinearSolvers> // pour le gradient conjugué

Vector Solve (const Matrix &A, const Vector &b, TYPE type)
{
    int size = mesh->GetNumberOfTotalPoints ();
    Vector sol (size);

    if (type == EXPLICIT) // U = Ab : produit matrice vecteur avec U = sol
    {
        sol = A * b;
    }
    if (type == IMPLICIT) // AU = b : gradient conjugé
    {
        ConjugateGradient<Matrix, Lower|Upper> cg; // déclare la fonction gradient conjugé
        cg.compute(A); // Calcule la matrice A
        sol = cg.solve(b) // Solution du problème
        sol = cg.solve(b) // met à jour la solution finale
    }

    return sol;
}
