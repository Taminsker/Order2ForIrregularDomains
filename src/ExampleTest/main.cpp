#include "../O2FID/O2FID.h"

#include <iostream>

double phi (Point a); // fonction levelset
double f (Point a); // fonction de second membre
double g (Point a); // fonction des conditions aux limites

double f_chaleur (Point a, double t); // fonction de second membre equation chaleur
double g_chaleur (Point a, double t); // fonction des conditions aux limites equation chaleur

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    //!!
    /// -------- Problème de Laplace  -------
    //!!
    {
        // ---- Création du maillage ---
        Mesh maillage = Mesh (); // Création d'un maillage vide

        maillage.SetBounds (Point(), Point (1, 1, 1)); // Configuration des points définissant le maillage ie [0, 1] ^3
        maillage.Set_Nx (3); // Le maillage aura 3 points dans la direction x
        maillage.Set_Ny (3); // 3 points dans la direction y
        maillage.Set_Nz (4); // et 4 points dans la direction z

        //! OBLIGATOIRE
        maillage.Build ();// Génération du maillage selon les paramètres précécdement choisis

        /// Optionnel
        //mesh.Print (); // Affichage des points
        ///

        // ---- Frontiere de Omega ---
        Border frontiere = Border (&maillage); // Création de la frontière
        frontiere.SetLevelSet (&phi); // Configuration de la fonction levelset utilisée
        std::vector <int> liste_index_points_frontiere = frontiere.MakeListOfIndexPoints (); // Création de la liste d'index de points sur la frontière ie les points de zéros de phi

        // ---- Construction de la matrice ---
        MatrixBuilder construcMatrice = MatrixBuilder (&maillage); // Création du monsieur qui va construire la matrice
        construcMatrice.Set3D (); // On lui précise que la matrice sera un problème 3D
        Matrix matrice = construcMatrice.MakeLaplaceEquation (); // On construit la matrice du problème de Laplace

        // ---- Création du Terme source ---
        FunToVec construTermeSource = FunToVec (&maillage);
        Vector secondMembre = construTermeSource.Make (&f); // On construit le second membre associé à la fonction f

        // ---- Imposition des conditions de Dirichlet ---
        Impose imposeur = Impose (&maillage); // Création du monsieur qui va imposer les conditions aux bords
        imposeur.SetIndexPoints (liste_index_points_frontiere); // On lui dit sur quels noeuds il va devoir imposer les conditions
        secondMembre -= imposeur.MakeDirichlet (&matrice, &g); // On soustrait on second membre les conditions qu'on impose

        // ---- Résolution du système en implicite ---
        Solver resolveur = Solver (&maillage); // Création du monsieur qui va résoudre
        resolveur.SetSolverToImplicit (); // On lui dit que ça va être une résolution en implicite
        Vector solution_numerique = resolveur.Solve (matrice, secondMembre); // On résout en précisant la matrice et le second membre

        // ---- Création du vecteur de solution analytique ---
        FunToVec construSolAna = FunToVec (&maillage);
        Vector solution_analytique = construSolAna.Make (&f); // On construit le second membre associé à la fonction f

        // ---- Calcul des erreurs ---
        ErrorsBuilder construcErreurs = ErrorsBuilder (&maillage); // On construit le monsieur qui va calculer toutes nos erreurs
        construcErreurs.SetVectorNumerical (&solution_numerique); // On lui donne le vecteur de la solution numérique
        construcErreurs.SetVectorAnalytical (&solution_analytique); // On lui donne le vecteur de la solution analytique
        double erreurL2 = construcErreurs.GetErrorL2 (); // On récupère l'erreur L2
        double erreurLinf = construcErreurs.GetErrorLinf (); // On récupère l'erreur Linf
        double erreurRela = construcErreurs.GetErrorRela (); // On récupère l'erreur relative
        Vector erreurAbs = construcErreurs.GetErrorAbs (); // On récupère le vecteur d'erreurs en valeur absolue

        // --- Écriture des fichiers ---
        Writer ecrivain = Writer (&maillage); // Création du monsieur qui va sauver les fichiers
        ecrivain.SetFilename ("Laplace_3D"); // On lui dit sous quel nom le sauver
        ecrivain.WriteBothDomainsOff (); // On lui dit d'écrire que ce qui est dans le domaine Omega 1
        ecrivain.SetCurrentIteration (0); // Il n'y a pas d'itération sur le temps donc 0
        ecrivain.SetVectorNumerical (&solution_numerique); // On lui donne la solution numerique
        ecrivain.SetVectorAnalytical (&solution_analytique); // On lui donne la solution analytique
        ecrivain.SetVectorErrorAbs (&erreurAbs); // On lui donne le vecteur des erreurs absolues
        ecrivain.WriteNow (); // Maintenant on lui demande d'écrire les fichiers !!

        // --- Affichage dans le terminal
        std::cout << "O2FID : Équation de la Laplace" << std::endl;
        std::cout << "\tDim = 3D" << std::endl;
        std::cout << "\tOuvert = [" << maillage.GetBounds () [0] << "]x[" << maillage.GetBounds ()[1] << "]" << std::endl;
        std::cout << "\tNx = " << maillage.Get_Nx () << "; hx = " << maillage.Get_hx () << std::endl;
        std::cout << "\tNy = " << maillage.Get_Ny () << "; hy = " << maillage.Get_hy () << std::endl;
        std::cout << "\tNz = " << maillage.Get_Nz () << "; hz = " << maillage.Get_hz () << std::endl;
        std::cout << "\tNombre total de points = " << maillage.GetNumberOfTotalPoints () << std::endl;
        std::cout << "\n\tMoyenne erreur abs = " << erreurAbs.mean () << std::endl;
        std::cout << "\tErreur rela = " << erreurRela << std::endl;
        std::cout << "\tErreur L2 = " << erreurL2 << std::endl;
        std::cout << "\tErreur Linf = " << erreurLinf << std::endl;
        std::cout << "\nFichiers d'écriture : Laplace_3D" << std::endl;
    }


// ------------------------------------------------------------------------------------------------


    //!!
    /// -------- Équation de la chaleur  -------
    //!!
    {
        double coeffChaleur = 14.5; // Exemple
        double dt = 1e-6; // Exemple
        double T = 3; // Exemple temps limite
        int nombreIteration = int (T / dt); // Nombre d'itération temporelle


        // ---- Création du maillage ---
        Mesh maillage = Mesh ();

        maillage.SetBounds (Point(), Point (1)); // [0, 1]

        //! OBLIGATOIRE
        maillage.Build ();

        // ---- Frontiere de Omega ---
        Border frontiere = Border (&maillage);
        frontiere.SetLevelSet (&phi);
        std::vector <int> liste_index_points_frontiere = frontiere.MakeListOfIndexPoint ();

        // ---- Construction de la matrice ---
        MatrixBuilder construcMatrice = MatrixBuilder (&maillage);
        construcMatrice.Set1D (); // On lui précise que la matrice sera un problème 1D
        construcMatrice.SetDeltaT (dt); // On lui précise le dt
        Matrix matrice = construcMatrice.MakeHeatEquation (coeffChaleur); // On construit la matrice du problème de la chaleur

        // ---- Création du Terme source ---
        FunToVec construTermeSource = FunToVec (&maillage);

        // ---- Imposition des conditions de Neumann ---
        Impose imposeur = Impose (&maillage);
        imposeur.SetIndexPoints (liste_index_points_frontiere);

        // ---- Résolution du système en explicite ---
        Solver resolveur = Solver (&maillage);
        resolveur.SetSolverToExplicit (); // On lui dit que ça va être une résolution en explicite

        // ---- Création du vecteur de solution analytique ---
        FunToVec construSolAna = FunToVec (&maillage);

        // ---- Calcul des erreurs
        ErrorsBuilder construcErreurs = ErrorsBuilder (&maillage);

        // --- Écriture des fichiers ---
        Writer ecrivain = Writer (&maillage);
        ecrivain.SetFilename ("Chaleur_1D"); // On lui dit sous quel nom le sauver


        // --- Affichage dans le terminal
        std::cout << "O2FID : Équation de la Laplace" << std::endl;
        std::cout << "\tDim = 3D" << std::endl;
        std::cout << "\tOuvert = [" << maillage.GetBounds () [0] << "]x[" << maillage.GetBounds ()[1] << "]" << std::endl;
        std::cout << "\tNx = " << maillage.Get_Nx () << "; hx = " << maillage.Get_hx () << std::endl;
        std::cout << "\tNy = " << maillage.Get_Ny () << "; hy = " << maillage.Get_hy () << std::endl;
        std::cout << "\tNz = " << maillage.Get_Nz () << "; hz = " << maillage.Get_hz () << std::endl;
        std::cout << "\tNombre total de points = " << maillage.GetNumberOfTotalPoints () << std::endl;

        // ---- Résolution ---
        for (int i = 0; i < nombreIteration ; i++)
        {
            std::cout << "\nOn résout pour le temps t = " << double (i) * dt << std::endl;

            // ---- Création du Terme source ---
            Vector secondMembre = construTermeSource.Make (&f);

            // ---- Imposition des conditions de Neumann ---
            secondMembre -= imposeur.MakeNeumann (&matrice, &g); // On soustrait on second membre les conditions qu'on impose

            // ---- Résolution du système en explicite ---
            Vector solution_numerique = resolveur.Solve (matrice, secondMembre);

            // ---- Création du vecteur de solution analytique ---
            Vector solution_analytique = construSolAna.Make (&f);

            // ---- Calcul des erreurs
            construcErreurs.SetVectorNumerical (&solution_numerique);
            construcErreurs.SetVectorAnalytical (&solution_analytique);
            double erreurL2 = construcErreurs.GetErrorL2 ();
            double erreurLinf = construcErreurs.GetErrorLinf ();
            double erreurRela = construcErreurs.GetErrorRela ();
            Vector erreurAbs = construcErreurs.GetErrorAbs ();

            // --- Écriture des fichiers ---
            ecrivain.SetCurrentIteration (i); // Il n'y a pas d'itération sur le temps donc 0
            ecrivain.SetVectorNumerical (&solution_numerique);
            ecrivain.SetVectorAnalytical (&solution_analytique);
            ecrivain.SetVectorErrorAbs (&erreurAbs);
            ecrivain.WriteNow ();

            std::cout << "\tMoyenne erreur abs = " << erreurAbs.mean () << std::endl;
            std::cout << "\tErreur rela = " << erreurRela << std::endl;
            std::cout << "\tErreur L2 = " << erreurL2 << std::endl;
            std::cout << "\tErreur Linf = " << erreurLinf << std::endl;
        }
    }

    return 0;
}

double phi (Point a)
{
    return 0. * (a.x + a.y + a.z);
}

double f (Point a)
{
    return 0. * (a.x + a.y + a.z);
}

double g (Point a)
{
    return 0. * (a.x + a.y + a.z);
}

double f_chaleur (Point a, double t)
{
    return 0. * (a.x + a.y + a.z) + t;
}

double g_chaleur (Point a, double t)
{
    return 0. * (a.x + a.y + a.z) + t;
}
