/** @file toolbox.cpp */

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

    vec->conservativeResize (N);

    for (int i = G; i < N; ++i)
    {
        double sum = 0.;
        std::vector<Point*> neigh = mesh->GetPoint(i)->GetListNeighbours ();

        for (Point* p_n : neigh)
        {
            sum += vec->coeff (p_n->GetGlobalIndex ());
        }

        vec->coeffRef (i) = sum / double(neigh.size ());
    }

    return;
}

void Extrapole (Mesh* mesh, std::vector<Point*>* vec)
{
    size_t N = size_t (mesh->GetNumberOfTotalPoints ());
    size_t G = size_t(mesh->GetNumberOfCartesianPoints ());

    vec->resize (N);
//    std::vector<Point*> R(N);

//    for(size_t i = 0; i < G; ++i)
//        R.at (i) = new Point (*vec->at (i));

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

        if (vec->at (i) != nullptr)
            delete vec->at (i);

        vec->at (i) = new Point(*sum / double(neigh.size () - size_t(count)));

        delete sum;

    }
    return;
}

void ExtendToNeighbours(Mesh* mesh, std::vector<int>* vec)
{
    auto R = *vec;
    size_t N = vec->size ();

    for (size_t i = 0; i < N; ++i)
    {
        auto neigh1 = mesh->GetPoint (vec->at (i))->GetListNeighbours ();

        for (auto p : neigh1)
        {
            auto neigh2 = p->GetListNeighbours ();

            for (auto n : neigh2)
            {
                auto neigh3 = n->GetListNeighbours ();

                for (auto v : neigh3)
                    R.push_back (v->GetGlobalIndex ());

                R.push_back (n->GetGlobalIndex ());
            }

            R.push_back (p->GetGlobalIndex ());
        }
    }

    std::sort(R.begin(), R.end());
    R.erase(std::unique(R.begin(), R.end()), R.end());

    *vec = R;

    return;
}


