#include "impose.h"

Vector ImposeDirichlet (Mesh * mesh,
                        Matrix * A,
                        double (*g) (Point, double),
                        std::vector <int> listIndex,
                        double t,
                        std::vector <double> secondMember)
{
  for (int i : listIndex) {A.row (i) *= 0.;} // On met la i-ème ligne à 0

  A = A.transpose (); // Pour pouvoir manipuler les colonnes de A

  for (int i : listIndex) {
    // On déplace au 2nd membre les apparitions de P_i avec valeur imposée g(P_i)
    secondMember -= g(*(mesh->GetPoint (i)), t) * A.row (i).transpose ();

    A.row (i) *= 0.;
    A.coeffRef (i,i) = 1.;

    secondMember [i] = g(*(mesh->GetPoint (i)), t);
  }

  A = A.transpose ();
  A.pruned ();
}


Vector ImposeNeumann (Mesh *mesh,
                      Matrix * sparsematrix,
                      double (*g) (Point, double),
                      double (*phi) (Point, double),
                      std::vector <int> listIndex,
                      double t)
{

    // Récupère le nombre points sur le bord
    size_t NBorder = listIndex.size ();

    // Créer le vecteur de conditions aux bords
    Vector cond_border (mesh->GetNumberOfTotalPoints ());

    // Boucle sur les indices des points sur le bords
    for (size_t i = 0; i < NBorder; i++)
    {
        // Récupère l'indice du point courant
        int index_p = listIndex.at (i);

        // Récupère le point p
        Point * p = mesh->GetPoint (index_p);

        // Récupère la liste des voisins du point courant
        std::vector <Point*> neigh = p->GetListNeighbours ();

        // Créer un vecteur des phi des voisins
        Vector phi_neighbour (neigh.size ());
        for (size_t j = 0; j < neigh.size (); ++j)
            phi_neighbour (int (j)) = phi (*neigh.at (j), t);

        // Stocke l'index du point qui est dans la zone (+) du domaine
        int index_n;
        // Récupère l'indice du point de la normale
        double coeff = phi_neighbour.maxCoeff (&index_n); // Récupère l'indice du phi max

        if (coeff <= 0) // Oups il semblerait qu'aucun des points de voisinages ne soit à l'exterieur du domaine
        {
            std::cout << "Erreur le point : \n" << i << p
                      << " du bord semble poser problème ! On le skip...\n\n";
            continue;
        }

        // Récupère le point "normal"
        int indexGlobal_p_normal = neigh.at (size_t (index_n))->GetGlobalIndex ();
        Point * p_normal = mesh->GetPoint (indexGlobal_p_normal);
        Point diff = (*p_normal - *p);

        if (TYPE_INTERPOLATION == DEGRE_1)
        {
            // u_p = u_N - g (x_p) * h
            double h = std::max (diff.x, std::max (diff.y, diff.z)); // Définition de la distance entre le point p et p_normal
            double g_p = g (*p, t); // Définition de la valeur g (x_p)

            sparsematrix->row (index_p) *= 0.; // Met la ligne_p à 0
            cond_border += g_p * h * sparsematrix->col (index_p); // Passe les conditions aux limites de l'autre coté

            sparsematrix->col (index_n) += sparsematrix->col (index_p); // Ajoute colonne_p à la colonne_n (correspond à remplace u_p par u_n)

            sparsematrix->col (index_p) *= 0.; // Met la colonne_p à 0
            sparsematrix->coeffRef (index_p, index_p) = 1.; // Ajoute le coefficient diagonal à 1;
        }
    }

    sparsematrix->pruned ();
    sparsematrix->makeCompressed ();
    return cond_border;
}
