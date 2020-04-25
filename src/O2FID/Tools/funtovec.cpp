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
        Point p_i = *mesh->GetPoint (i);

        /* On affecte la valeur de f à ce point. */
        Fvec.coeffRef (i) = f(p_i, t);
    }

    return Fvec;
}

Vector FunToVec (Mesh * mesh, double value)
{
    //    /* Nombre de point total dans le mesh à calculer. */
    int n = mesh->GetNumberOfTotalPoints ();
    //    /* Création d'un vecteur de n points. */
    Vector Fvec (n);

    //    for (int i = 0; i < n; i++)
    //    {
    //        Fvec (i) = value;
    //    }

    return Fvec.setConstant (value);
}

Vector FunToVec (Mesh * mesh,
                 Point (*f) (double, double),
                 Point (*fprim) (double, double),
                 size_t numOfPtsBorder,
                 double t,
                 double a,
                 double b)
{
    std::cout << "# FunToVec parametric with " << numOfPtsBorder << " on border." << std::endl;

    int numPoints = mesh->GetNumberOfTotalPoints ();
    // Vecteur de sortie
    Vector out (numPoints);
    out.setZero ();

    if (numOfPtsBorder <= 0)
        return out;

    double h = double (b - a) / double (numOfPtsBorder);

    std::vector<Point *> position (numOfPtsBorder);

    for (size_t i = 0; i < numOfPtsBorder; ++i)
    {
        position.at (i) = new Point (f (a + i * h, t));
        position.at (i)->SetGlobalIndex (int (i));
    }

    for (int idx = 0; idx < numPoints; ++idx)
    {
        Point * p = mesh->GetPoint (idx);

        Point * p_min = position.at (0);
        double dist = EuclidianDist (*p, *p_min);

        for (Point * p_curve : position)
        {
            double d = EuclidianDist (*p, *p_curve);

            if (d < dist)
            {
                p_min = p_curve;
                dist = d;
            }
        }

        Point p_tangent = fprim (a + p_min->GetGlobalIndex () * h, t);
        Point normal_out = {-p_tangent.y, p_tangent.x, 0.};

        out.coeffRef (idx) = (normal_out | (*p_min - *p));

        std::cout << "\r" << INDENT << "The minimum was found for point " << idx+1 << "/" << numPoints << " at edge point " << p_min->GetGlobalIndex () << ".          " << std::flush;
    }

    std::cout << "\r" << INDENT << "The minimum was found for all the points.                                 " << std::endl;

    // Suppression du vecteur de position
    for (Point * p : position)
        delete p;

    position.clear ();

    // Fin de l'algo
    std::cout << std::endl;

    return out;
}
