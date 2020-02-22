/** \file solver.h */

#ifndef SOLVER_H
#define SOLVER_H

#include "../Base/abstractbase.h"
#include "../Data/data.h"
#include "../Mesh/mesh.h"

/**
 * @brief La classe Solver implémente un solver implicit et un solveur explicite
 */
class Solver : public AbstractBase
{
public:
    typedef enum {
        EXPLICIT,
        IMPLICIT
    } TYPE;

    /**
     * @brief Conctructeur
     */
    Solver (Mesh * mesh);

    /**
     * @brief Destructeur
     */
    ~Solver ();

    /**
     * @brief Configure le solveur sur un solveur implicite.
     */
    void SetSolverToImplicit ();

    /**
     * @brief Configure le solveur sur un solveur explicite.
     */
    void SetSolverToExplicit ();

    /**
     * @brief Solve
     * @param A Matrix *
     * @param b Vector *
     * @return
     */
    Vector Solve (const Matrix &A, const Vector &b);

protected:
    int m_type = EXPLICIT;
};

#endif // SOLVER_H
