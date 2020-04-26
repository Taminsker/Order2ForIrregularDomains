
#include "../Data/data.h"
#include "../Mesh/mesh.h"
#include "../Data/datatypedefinitions.h"
#include "toolbox.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

int Remainder (int dividend, int divisor)
{
    while (dividend >= divisor || dividend < 0)
        dividend += (dividend < 0 ? 1 : -1) * divisor;

    return dividend;
}

int Quotient (int dividend, int divisor)
{
    int quotient = 0;
    while (dividend >= divisor || dividend < 0)
    {
        int a = (dividend < 0 ? 1 : -1);
        quotient += a;
        dividend += a * divisor;
    }

    return quotient;
}

std::vector<double> Order (std::vector<double> err, std::vector<double> h)
{
    size_t N = std::min (err.size (), h.size ());

    std::vector<double> order (N, 0.);

    for (size_t i = 1; i < N; ++i)
        order.at (i) = (std::log (err.at (i)) - std::log( err.at (i-1))) / (std::log (h.at (i)) - std::log(h.at (i-1)));

    return order;
}

void Extrapole (Mesh* mesh, Vector* vec)
{
    int N = mesh->GetNumberOfTotalPoints ();
    int G = mesh->GetNumberOfCartesianPoints ();
    Vector R(N);
    R.setZero ();

    for(int i = 0; i < G; ++i)
        R.coeffRef (i) = vec->coeffRef (i);

    for (int i = G; i < N; ++i)
    {
        double sum = 0.;
        std::vector<Point*> neigh = mesh->GetPoint(i)->GetListNeighbours ();

        for (Point* p_n : neigh)
        {
            sum += vec->coeffRef (p_n->GetGlobalIndex ());
        }

        R.coeffRef (i) = sum / double(neigh.size ());
    }

    *vec = R;

    return;
}

void Extrapole (Mesh* mesh, std::vector<Point*>* vec)
{
    size_t N = size_t (mesh->GetNumberOfTotalPoints ());
    size_t G = size_t(mesh->GetNumberOfCartesianPoints ());

    std::vector<Point*> R(N);

    for(size_t i = 0; i < G; ++i)
        R.at (i) = new Point (*vec->at (i));

    for (size_t i = G; i < N; ++i)
    {
        Point* sum = new Point();

        std::vector<Point*> neigh = mesh->GetPoint(int(i))->GetListNeighbours ();

        int count = 0;
        for (Point* p_n : neigh)
        {
            if (p_n->GetGlobalIndex () < int(vec->size ()))
                *sum = *sum + *vec->at (size_t (p_n->GetGlobalIndex ()));
            else
                count++;
        }

        R.at (i) = new Point(*sum / double(neigh.size () - size_t(count)));
    }

    AutoClearVector(vec);

    *vec = R;

    return;
}

