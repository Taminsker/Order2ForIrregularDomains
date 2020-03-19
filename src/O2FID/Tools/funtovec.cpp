#include "funtovec.h"

Vector FunToVec (Mesh * mesh, double (*f) (Point, double), double t)
{
    /* Nombre de point total dans le mesh à calculer. */
    int n = mesh->GetNumberOfTotalPoints ();
    /* Création d'un vecteur de n points. */
    Vector Fvec (n);

    for (int i = 0; i < n; i++)
    {
        /* On récupère le i-ème point du maillages. */
        Point p_i = mesh->operator() (i);

        /* On affecte la valeur de f à ce point. */
        Fvec (i) = f(p_i, t);
    }

    return Fvec;
}

Vector FunToVec (Mesh * mesh, double value)
{
//    /* Nombre de point total dans le mesh à calculer. */
//    int n = mesh->GetNumberOfTotalPoints ();
//    /* Création d'un vecteur de n points. */
//    Vector Fvec (n);

//    for (int i = 0; i < n; i++)
//    {
//        Fvec (i) = value;
//    }

//    return Fvec;

    int n = mesh->GetNumberOfTotalPoints ();
    Vector Fvec (n);
    Fvec.setConstant (value);
    return Fvec;
}
