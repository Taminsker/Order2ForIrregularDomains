#include "border.h"

#include <stdio.h>
#include <stdlib.h>


std::vector <int> MakeListOfIndexPoints (Mesh * mesh, double (*phi) (Point, double), double t)
{
  typedef Point<double> Point;

  double eps = 1e-10; // Seuil sous lequel on considère les valeurs nulles

  std::vector<int> ListOfIndexPoints;

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; i < Ny; j++) {
      for (int k = 0; k < Nz; k++) {

        // On note P_indice le point suivant dans la direction représentée par l'indice
        Point p = Point (i, j, k);
        Point p_i = Point (i + 1, j, k); // Le point en i+1
        Point p_j = Point (i, j + 1, k); // Le point en j+1
        Point p_k = Point (i, j, k + 1); // Le point en k+1

        // On note val_indice la valeur de phi au point P_indice
        val = phi (p, t);
        val_i = phi (p_i, t); // La valeur en i+1
        val_j = phi (p_j, t); // La valeur en j+1
        val_k = phi (p_k, t); // La valeur en k+1

        // Si le noeud (i,j,k) du maillage est déjà lui-même un point du bord
        if (fabs(val) < eps)
        {
          p.SetLocate (ON_BORDER_OMEGA);
        }

        // Si le bord coupe une arête de direction (Ox), on définit un point intermédiaire
        else if ( (i + 1) < Nx && fabs(val_i) >= eps && (val * val_i) < 0 )
        {
          double di = fabs(val) / ( fabs(val) + fabs(val_i) );
          int num_global = mesh->AddPointOnBorder ( p + Point (di * hx, 0., 0.) );
        }

        // Si le bord coupe une arête de direction (Oy), on définit un point intermédiaire
        else if ( (j + 1) < Nx && fabs(val_j) >= eps && (val * val_j) < 0 )
        {
          double dj = fabs(val) / ( fabs(val) + fabs(val_j) );
          int num_global = mesh->AddPointOnBorder ( p + Point (0., dj * hy, 0.) );
        }

        // Si le bord coupe une arête de direction (Oz), on définit un point intermédiaire
        else if ( (k + 1) < Nx && fabs(val_k) >= eps && (val * val_k) < 0 )
        {
          double dk = fabs(val) / ( fabs(val) + fabs(val_k) );
          int num_global = mesh->AddPointOnBorder ( p + Point (0., 0., dk * hz) );
        }

      }
    }
  }

  return mesh->GetListOfIndexPoints ();

}
