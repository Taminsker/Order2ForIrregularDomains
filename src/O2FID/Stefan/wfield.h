#ifndef WFIELD_H
#define WFIELD_H

#include "../Mesh/mesh.h"
#include "../Data/data.h"
#include "../Data/datatypedefinitions.h"

#include "../Tools/matrix.h"
#include "../Tools/impose.h"
#include "../Tools/solver.h"
#include "../Tools/border.h"

#include "../Toolbox/differencefinite.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>



std::vector<Point*> GetWField(Mesh* mesh, Vector* phi, Vector* sol, std::vector<int>* idxsBorder, double h0 = 1);
std::vector<Point*> BuildWOnBorder (Mesh* mesh, Vector* phi, Vector* sol, std::vector<int>* idxsBorder, double h0);
void ExtendWToAllDomain (Mesh* mesh, std::vector<Point*>* W, std::vector<int>* idxsBorder);

#endif // WFIELD_H
