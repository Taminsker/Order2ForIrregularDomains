#ifndef INTERFACE_H
#define INTERFACE_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"
#include "../Data/datatypedefinitions.h"

#include "../Tools/matrix.h"
#include "../Tools/impose.h"
#include "../Tools/solver.h"
#include "../Tools/border.h"

#include "../Toolbox/differencefinite.h"
#include "../Toolbox/toolbox.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>


Vector IteratePhi (Mesh* mesh, std::vector<Point*>* W, double dt, Vector* phi_n, Vector* phi_n_1 = nullptr, Vector* phi_n_2 = nullptr);

Vector Weno(Mesh* mesh, Vector* u, std::vector<Point*>* W);

void ReInitPhi (Mesh* mesh, Vector* phi, std::vector<int>* idxsBorder, double dt);


#endif // INTERFACE_H
