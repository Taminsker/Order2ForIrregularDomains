/** @file border.cpp */

#include "border.h"


std::vector <int> MakeBorderPoints (Mesh * mesh, Vector * phi_list)
{
    std::cout << "# MakeBorderPoints" << std::endl;

    double eps = 1e-15; // Seuil sous lequel on considère les valeurs nulles

    int numberOfPointsOnGrid = mesh->GetNumberOfTotalPoints ();

    int GridPoints = 0;
    int NewPoint = 0;

    double maximum = std::max (mesh->Get_hx (), std::max (mesh->Get_hy (), mesh->Get_hz ())) + 1e-3;

    for (int i = 0; i < numberOfPointsOnGrid; i++)
    {

        Point * p_curr = mesh->GetPoint (i);

        double phi_curr = phi_list->coeff (i);

        if (fabs(phi_curr) < eps)
        {
            // On tag le point
            p_curr->SetLocate (ON_BORDER_OMEGA);
            GridPoints++;
            continue;
        } else
            p_curr->SetLocate ((phi_curr < 0. ? ON_DOMAIN_INTERN_OMEGA : ON_DOMAIN_EXTERN_OMEGA));

        for (Point * p_neigh : p_curr->GetListNeighbours ())
        {
            // Le point p_neigh est un point rajouté
            if (p_neigh->GetGlobalIndex () > mesh->GetNumberOfCartesianPoints ())
                continue;

            double phi_neigh = phi_list->coeff (p_neigh->GetGlobalIndex ());

            // ordre croissant des indices
            if (p_neigh->GetGlobalIndex () >= p_curr->GetGlobalIndex ())
                continue;

            // p_neigh est un point de bord
            if (fabs (phi_neigh) < eps)
                continue;

            // Même signe
            if (phi_curr * phi_neigh > 0.)
                continue;

            // trop éloigné /!\ attention symétrie
            if (EuclidianDist (*p_curr, *p_neigh) > maximum)
                continue;

            double dist = fabs(phi_curr) / (fabs(phi_curr) + fabs(phi_neigh));

            Point * p_new = new Point (*p_curr + dist * (*p_neigh - *p_curr));

            // Suppression des voisins
            p_curr->RemoveThisNeighbourPoint (p_neigh);
            p_neigh->RemoveThisNeighbourPoint (p_curr);

            // Ajout des voisins
            p_curr->AddPointNeighbour (p_new);
            p_neigh->AddPointNeighbour (p_new);

            // Ajout des voisins
            p_new->AddPointNeighbour (p_curr);
            p_new->AddPointNeighbour (p_neigh);

            // Ajout du point au maillage
            mesh->AddPointOnBorder (p_new);

            // Vertex
            mesh->MakeACellFromListPoints ({p_new});

            NewPoint++;
        }
    }

    std::cout << INDENT << GridPoints << " points of the cartesian grid are considered on the edge of the domain." << std::endl;
    std::cout << INDENT << NewPoint << " points have been added on edges." << std::endl;

    std::cout << std::endl;

    return mesh->GetListOfIndexPoints ();
}

