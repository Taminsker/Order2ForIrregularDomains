#ifndef CHALEUR1D
#define CHALEUR1D

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

Eigen::VectorXd mesh (double a, double b, int n);

Eigen::SparseMatrix<double> matrix (double a, double b, int n, double dt,
                                    double D);

#endif
