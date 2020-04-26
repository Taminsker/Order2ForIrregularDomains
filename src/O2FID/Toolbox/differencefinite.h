#ifndef DIFFFINITE_H
#define DIFFFINITE_H

#include <vector>
#include<algorithm>

#include "toolbox.h"

struct DFIdxCoeff
{
    std::vector<int> idxs = {};
    std::vector<double> coeffs = {};

};

struct DFOrders
{
    DFIdxCoeff Order1;
    DFIdxCoeff Order2;
    DFIdxCoeff Order3;
    DFIdxCoeff Order4;
    DFIdxCoeff Order5;
    DFIdxCoeff Order6;
    DFIdxCoeff Order7;
    DFIdxCoeff Order8;
};

struct DerivativeStruct
{
    DFOrders Backward;
    DFOrders Central;
    DFOrders Forward;
};

struct DFStruct
{
    DerivativeStruct Derivative_1;
    DerivativeStruct Derivative_2;

    DFStruct();
};

DFOrders Reverse (DFOrders O);
DFIdxCoeff Reverse (DFIdxCoeff D);
#endif // DIFFFINITE_H
