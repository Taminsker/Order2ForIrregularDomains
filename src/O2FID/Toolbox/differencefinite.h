#ifndef DIFFFINITE_H
#define DIFFFINITE_H

#include <vector>
#include<algorithm>

#include "toolbox.h"

typedef enum
{
    ORDER_1_FORWARD,
    ORDER_2_FORWARD,
    ORDER_3_FORWARD,
    ORDER_4_FORWARD,
    ORDER_5_FORWARD,
    ORDER_6_FORWARD,
    ORDER_7_FORWARD,
    ORDER_8_FORWARD,
    ORDER_1_CENTRAL,
    ORDER_2_CENTRAL,
    ORDER_3_CENTRAL,
    ORDER_4_CENTRAL,
    ORDER_5_CENTRAL,
    ORDER_6_CENTRAL,
    ORDER_7_CENTRAL,
    ORDER_8_CENTRAL,
    ORDER_1_BACKWARD,
    ORDER_2_BACKWARD,
    ORDER_3_BACKWARD,
    ORDER_4_BACKWARD,
    ORDER_5_BACKWARD,
    ORDER_6_BACKWARD,
    ORDER_7_BACKWARD,
    ORDER_8_BACKWARD,
} ORDERS;


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

void DFOrderBuild (int der, ORDERS order, std::vector<int>* idxs, std::vector<double>* coeffs);

DFOrders Reverse (DFOrders O);
DFIdxCoeff Reverse (DFIdxCoeff D);
#endif // DIFFFINITE_H
