/** \file typedefinitions.h **/

#ifndef TYPEDEFINITIONS_H
#define TYPEDEFINITIONS_H

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <utility>

/**
 * @brief Abréviation du type Matrix
 */
typedef Eigen::SparseMatrix<double> Matrix;

/**
 * @brief Abréviation du type Vector
 */
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

#endif // TYPEDEFINITIONS_H
