#include "impose.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

Matrix BuildMatrixLaplaceEquation (Mesh * mesh)
{
/*
// Nombre de point totale
int n = mesh->GetNumberOfTotalPoints ();

// Pas selon les différentes directions spatiale en concidérant si 1D, 2D ou 3D
double hx = (mesh->Get_hx () < 1e-8) ? 1 : mesh->Get_hx ();
double hy = (mesh->Get_hy () < 1e-8) ? 1 : mesh->Get_hy ();
double hz = (mesh->Get_hz () < 1e-8) ? 1 : mesh->Get_hz ();
*/
}


Matrix BuildMatrixHeatEquation (Mesh * mesh, double dt, double D, int TYPE)
{
//  // On récupère le nombre de points et les pas d'espace

//  int Nx = mesh->Get_Nx ();
//  int Ny = mesh->Get_Nx ();
//  int Nz = mesh->Get_Nz ();

//  int N = mesh->GetNumberOfTotalPoints ();

//  double hx = mesh->Get_hx ();
//  double hy = mesh->Get_hx ();
//  double hz = mesh->Get_hz ();

//  double a = 0., b = 0., c = 0., d = 0.;

//  // On définit les coefficients pour implicite par défaut

//  if (hx < 1e-10)
//  {
//    b = - (D * dt) / (hx * hx);
//  }
//  if (hy < 1e-10)
//  {
//    c = - (D * dt) / (hy * hy);
//  }
//  if (hz < 1e-10)
//  {
//    d = - (D * dt) / (hz * hz);
//  }
//  a = 1. - 2.*b - 2.*c - 2.*d;

//  // On définit un alias pour les triplets

//  typedef Eigen::Triplet<double> Triplet;

//  // On définit la liste de triplets et on prépare la place

//  std::vector<Triplet> ListOfTriplets;
//  ListOfTriplets.reserve(7 * N);

//  // On remplit la liste de triplets

//  for (int i = 0; i < N; ++i) {

//    int surDiag1 = i + 1;
//    int surDiag2 = i + Nx;
//    int surDiag3 = i + (Nx * Ny);

//    // La diagonale principale des a

//    ListOfTriplets.push_back(Triplet(i,i,a));

//    // La sur-diagonale des b, décalée de 1 (si possible)

//    if (surDiag1 < N)
//    {
//      ListOfTriplets.push_back(Triplet(i,surDiag1,b));
//    }

//    // La sur-diagonale des c, décalée de Nx (si possible)

//    if (surDiag2 < N)
//    {
//      ListOfTriplets.push_back(Triplet(i,surDiag2,c));
//    }

//    // La sur-diagonale des d, décalée de Nx * Ny (si possible)

//    if (surDiag3 < N)
//    {
//      ListOfTriplets.push_back(Triplet(i,surDiag3,d));
//    }

//  }

//  // On construit la partie supérieure de la matrice à partir de la liste de triplets

//  Matrix A (N, N);
//  A.setFromTriplets(ListOfTriplets.begin(), ListOfTriplets.end());

//  // Si on résout en explicite, A devient -A + 2*identité

//  if (TYPE == EXPLICIT)
//  {
//    Matrix C (N, N);
//    C.setIdentity ();
//    A *= -1.;
//    A += 2. * C;
//  }

//  // On lui ajoute sa transposée pour avoir la matrice complète

//  A = A.selfadjointView<Upper> ();

//  return A;
}
