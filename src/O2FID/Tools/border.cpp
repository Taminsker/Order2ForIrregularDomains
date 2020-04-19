#include "border.h"


std::vector <int> MakeListOfIndexPoints (Mesh * mesh, Vector * phi_list, double t)
{
    double eps = 1e-10; // Seuil sous lequel on considère les valeurs nulles

    int numberOfPointsOnGrid = mesh->GetNumberOfTotalPoints ();

    for (int i = 0; i < numberOfPointsOnGrid; i++)
    {

        Point * p_curr = mesh->GetPoint (i);

        double phi_curr = phi_list->coeff (i);

        if (fabs(phi_curr) < eps)
        {
            // On tag le point
            p_curr->SetLocate (ON_BORDER_OMEGA);
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

            double dist = fabs(phi_curr) / (fabs(phi_curr) + fabs(phi_neigh));

            // Tests
//            double dist_p = std::sqrt((*p_curr - *p_neigh)| (*p_curr - *p_neigh));
//            bool b1 = (fabs(dist_p - mesh->Get_hx ()) > eps);
//            bool b2 = (fabs(dist_p - mesh->Get_hx ()) > eps);
//            bool b3 = (fabs(dist_p - mesh->Get_hx ()) > eps);

            Point * p_new = new Point (*p_curr + dist * (*p_neigh - *p_curr));

//            if (b1 && b2 && b3)
//            {
//                std::cout << "Add a point : " << std::endl;
//                std::cout << "\t p_curr (" << p_curr->GetGlobalIndex () << ")\t" << *p_curr << " phi : " << phi_curr << std::endl;
//                std::cout << "\t p_neigh (" << p_neigh->GetGlobalIndex () << ")\t" << *p_neigh << " phi : " << phi_neigh << std::endl;
//                std::cout << "\t p_new (" << mesh->GetNumberOfTotalPoints () << ")\t" << *p_new << std::endl;
//                std::cout << "\tdist(p_curr, p_neigh) = " << dist_p << std::endl;
//            }


            // Suppression des voisins
            p_curr->RemoveThisNeighbourPoint (p_neigh);
            p_neigh->RemoveThisNeighbourPoint (p_curr);

            // Ajout des voisins
            p_curr->AddPointNeighbour (p_new);
            p_neigh->AddPointNeighbour (p_new);

            // Ajout du point au maillage
            mesh->AddPointOnBorder (p_new);

            // Vertex
            mesh->MakeACellFromListPoints ({p_new});
        }
    }

    return mesh->GetListOfIndexPoints ();
}

