/** @file tools.h */
#ifndef TOOLS_H
#define TOOLS_H

/**
 * \defgroup Outils
 * @brief Ce module regroupe toutes les fonctions utiles pour le projet.
 *  * \code
 * #include "../O2FID/O2FID.h"
 * #include <iostream>
 * double phi (Point a, double t = 0); // fonction levelset
 * double f_laplace (Point a, double t = 0.); // fonction de second membre
 * double g_laplace (Point a, double t = 0.); // fonction des conditions aux limites
 *
 * double f_chaleur (Point a, double t); // fonction de second membre equation chaleur
 * double g_chaleur (Point a, double t); // fonction des conditions aux limites equation chaleur
 *
 * int main(int argc, char* argv[])
 * {
 *  (void)argc;
 *  (void)argv;
 *
 *  //!!
 *   /// -------- Problème de Laplace  -------
 *   //!!
 *   {
 *       // ---- Création du maillage ---
 *       Mesh maillage = Mesh (); // Création d'un maillage vide
 *
 *       maillage.SetBounds (Point(), Point (1, 1, 1)); // Configuration des points définissant le maillage ie [0, 1]^3
 *       maillage.Set_Nx (3); // Le maillage aura 3 points dans la direction x
 *       maillage.Set_Ny (3); // 3 points dans la direction y
 *       maillage.Set_Nz (4); // et 4 points dans la direction z
 *
 *      //! OBLIGATOIRE
 *      maillage.Build ();// Génération du maillage selon les paramètres précédemment choisis
 *
 *      /// Optionnel
 *      //mesh.Print (); // Affichage des points
 *      ///
 *
 *      // ---- Frontiere de Omega ---
 *       std::vector <int> liste_index_points_frontiere = MakeListOfIndexPoints (&maillage, &phi); // Création de la liste d'index de points sur la frontière ie les points de zéros de phi
 *
 *
 *      // ---- Construction de la matrice ---
 *       Matrix matrice = BuildMatrixLaplaceEquation (&maillage); // On construit la matrice du problème de Laplace
 *
 *      // ---- Création du Terme source ---
 *       Vector secondMembre = FunToVec (&maillage, &f_laplace); // On construit le second membre associé à la fonction f
 *
 *       // ---- Imposition des conditions de Dirichlet ----
 *      secondMembre -= ImposeDirichlet (&maillage, &matrice, &g_laplace, liste_index_points_frontiere); // On soustrait on second membre les conditions qu'on impose
 *
 *      // ---- Résolution du système en implicite ----
 *       Vector solution_numerique = Solve (matrice, secondMembre, IMPLICIT); // On résout en précisant la matrice et le second membre
 *
 *      // ---- Création du vecteur de solution analytique ----
 *       Vector solution_analytique = FunToVec (&maillage, &f_laplace); // On construit le second membre associé à la fonction f
 *
 *      // ---- Calcul des erreurs ----
 *
 *      double erreurL2 = GetErrorL2 (&maillage, solution_analytique, solution_numerique); // On récupère l'erreur L2
 *      double erreurLinf = GetErrorLinf (&maillage, solution_analytique, solution_numerique); // On récupère l'erreur Linf
 *      double erreurRela = GetErrorRela (&maillage, solution_analytique, solution_numerique); // On récupère l'erreur relative
 *      Vector erreurAbs = GetErrorAbs(&maillage, solution_analytique, solution_numerique); // On récupère le vecteur d'erreurs en valeur absolue
 *
 *      // ---- Écriture des fichiers ----
 *      Writer ecrivain = Writer (&maillage); // Création du monsieur qui va sauver les fichiers
 *      ecrivain.SetFilename ("Laplace_3D"); // On lui dit sous quel nom le sauver
 *      ecrivain.WriteBothDomainsOff (); // On lui dit d'écrire que ce qui est dans le domaine Omega 1
 *      ecrivain.SetCurrentIteration (0); // Il n'y a pas d'itération sur le temps donc 0
 *      ecrivain.SetVectorNumerical (&solution_numerique); // On lui donne la solution numerique
 *      ecrivain.SetVectorAnalytical (&solution_analytique); // On lui donne la solution analytique
 *      ecrivain.SetVectorErrorAbs (&erreurAbs); // On lui donne le vecteur des erreurs absolues
 *      ecrivain.WriteNow (); // Maintenant on lui demande d'écrire les fichiers !!
 *
 *     // ---- Affichage dans le terminal ----
 *      std::cout << "O2FID : Équation de la Laplace" << std::endl;
 *      std::cout << "\tDim = 3D" << std::endl;
 *      std::cout << "\tOuvert = [" << maillage.GetBounds () [0] << "]x[" << maillage.GetBounds ()[1] << "]" << std::endl;
 *      std::cout << "\tNx = " << maillage.Get_Nx () << "; hx = " << maillage.Get_hx () << std::endl;
 *      std::cout << "\tNy = " << maillage.Get_Ny () << "; hy = " << maillage.Get_hy () << std::endl;
 *      std::cout << "\tNz = " << maillage.Get_Nz () << "; hz = " << maillage.Get_hz () << std::endl;
 *      std::cout << "\tNombre total de points = " << maillage.GetNumberOfTotalPoints () << std::endl;
 *      std::cout << "\n\tMoyenne erreur abs = " << erreurAbs.mean () << std::endl;
 *      std::cout << "\tErreur rela = " << erreurRela << std::endl;
 *      std::cout << "\tErreur L2 = " << erreurL2 << std::endl;
 *      std::cout << "\tErreur Linf = " << erreurLinf << std::endl;
 *      std::cout << "\nFichiers d'écriture : Laplace_3D" << std::endl;
 *  }
 *
 *  return 0;
 * }
 *
 * double phi (Point a, double t)
 * {
 *   (void)t;
 *  return 0. * (a.x + a.y + a.z);
 * }
 *
 * double f_laplace (Point a, double t)
 * {
 *  (void)t;
 *  return 0. * (a.x + a.y + a.z);
 * }
 * \endcode
 */
#include "border.h"
#include "errors.h"
#include "funtovec.h"
#include "impose.h"
#include "matrix.h"
#include "solver.h"

#endif // TOOLS_H
