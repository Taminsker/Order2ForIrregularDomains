#include "tools.h"

Vector ImposeDirichlet (Mesh * mesh,
                        Matrix * sparsematrix,
                        double (*g) (Point, double),
                        std::vector <int> listIndex,
                        double t)
{

}

Vector ImposeNeumann (Mesh &mesh,
                      Matrix * sparsematrix,
                      double (*g) (Point, double),
                      double (*phi) (Point, double),
                      std::vector <int> listIndex,
                      double t)
{

    // Récupère le nombre points sur le bord
    int NBorder = listIndex.size ();

    // Créer le vecteur de conditions aux bords
    Vector cond_border (mesh.GetNumberOfTotalPoints ());

    // Boucle sur les indices des points sur le bords
    for (int i = 0; i < NBorder; i++)
    {
        // Récupère l'indice du point courant
        int index_p = listIndex.at (i);

        // Récupère le point p
        Point p = mesh (index_p);

        // Récupère la liste des voisins du point courant
        std::vector <Point> neigh = p.GetListNeighbours ();

        // Créer un vecteur des phi des voisins
        Vector phi_neigh (neigh.size ());
        for (int j = 0; j < neigh.size (); ++j)
            phi_neigh [j] = phi (neigh.at (j), t);

        // Stock l'index du point qui est dans la zone (+) du domaine
        int index_normal;
        // Récupère l'indice du point de la normale
        double coeff = phi_neigh.maxCoeff (&index_normal); // Récupère l'indice du phi max

        if (coeff <= 0) // Oups il semblerait qu'aucun des points de voisinages ne soit à l'exterieur du domaine
        {
            std::cout << "Erreur le point : \n" << i << p
                      << " du bord semble poser problème ! On le skip...\n\n";
            continue;
        }

        // Récupère le point "normal"
        Point p_normal = mesh (neigh.at (index_normal).GetGlobalIndex ());
        Point diff = (p_normal - p);

        // u_p = u_N - g (x_p) * h
        double h = std::max (diff.x, std::max (diff.y, diff.z)); // Définition de la distance entre le point p et p_normal
        double g_p = g (p, t); // Définition de la valeur g (x_p)

        sparsematrix->row (index_p) *= 0.; // Met la ligne_p à 0
        cond_border += g_p * h * sparsematrix->col (index_p); // Passe les conditions aux limites de l'autre coté

        sparsematrix->col (index_n) += sparsematrix->col (index_p); // Ajoute colonne_p à la colonne_n (correspond à remplace u_p par u_n)

        sparsematrix->col (index_p) *= 0.; // Met la colonne_p à 0
        sparsematrix->coeffRef (index_p, index_p) = 1.; // Ajoute le coefficient diagonal à 1;
    }

    sparsematrix->makeCompressed ();
    return;
}
