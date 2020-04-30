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

#include <vector>

class Field
{
public:
    std::vector<Point*> W;
    std::vector<Point*> Normals;
    std::vector<Point*> GradTemperature;
    std::vector<Point*> GradPhi;

    Field() :
        W(std::vector<Point*>()),
        Normals(std::vector<Point*>()),
        GradTemperature(std::vector<Point*>()),
        GradPhi(std::vector<Point*>())
    {}

    ~Field ()
    {
        AutoClearVector(&W);
        AutoClearVector(&Normals);
        AutoClearVector(&GradTemperature);
        AutoClearVector(&GradPhi);
    }
    Field& operator= (const Field& f)
    {
        AutoClearVector(&W);
        AutoClearVector(&Normals);
        AutoClearVector(&GradTemperature);
        AutoClearVector(&GradPhi);

        W = f.W;
        Normals = f.Normals;
        GradTemperature = f.GradTemperature;
        GradPhi = f.GradPhi;

        return *this;
    }
};

Field GetWField (Mesh* mesh, Vector* phi, Vector* sol, std::vector<int>* idxsBorder, double h0 = 1);

Field BuildWOnBorder (Mesh* mesh, Vector* phi, Vector* sol, std::vector<int>* idxsBorder, double h0);

void ExtendWToAllDomain (Field* field, Mesh* mesh, std::vector<int>* idxsBorder);

#endif // WFIELD_H
