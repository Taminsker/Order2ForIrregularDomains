#include "border.h"

#include <stdio.h>
#include <stdlib.h>


std::vector <int> MakeListOfIndexPoints (Mesh * mesh, double (*phi) (Point, double))
{
  typedef Point<double> Point;

  // On initialise un vecteur qui contiendra les points de la frontière du domaine
  int reserve = Nx * Ny;
  std::vector<Point> ListOfPoints;
  ListOfPoints.reserve(reserve);

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; i < Ny; j++) {

      // Si le noeud (i,j) du maillage est déjà lui-même un point du bord
      if ( phi(i * hx, j * hy) == 0 )
      {
        ListOfPoints.push_back(Point(i,j));
      }

      // Si le bord coupe une arête horizontale, on définit un point intermédiaire
      else if ( (i + 1) < Nx && phi(i * hx, j * hy) * phi((i + 1) * hx, j * hy) < 0 )
      {
        double di = fabs( phi(i * hx, j * hy) ) / ( fabs( phi(i * hx, j * hy) ) + fabs( phi((i + 1) * hx, j * hy) ) );
        ListOfPoints.push_back(Point( (i + di) * hx, j * hy ));
      }

      // Si le bord coupe une arête verticale, on définit un point intermédiaire
      else if ( (j + 1) < Ny && phi(i * hx, j * hy) * phi(i * hx, (j + 1) * hy) < 0 )
      {
        double dj = fabs( phi(i * hx, j * hy) ) / ( fabs( phi(i * hx, j * hy) ) + fabs( phi(i * hx, (j + 1) * hy) ) );
        ListOfPoints.push_back(Point( i * hx, (j + dj) * hy ));
      }

    }
  }

}
