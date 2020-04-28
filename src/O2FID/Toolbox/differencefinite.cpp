#include "differencefinite.h"


DFStruct::DFStruct ()
{
    // Derivée première
    // Forward
    Derivative_1.Forward.Order1.idxs = {0, 1};
    Derivative_1.Forward.Order1.coeffs = {-1., 1.};

    Derivative_1.Forward.Order2.idxs = {0, 1, 2};
    Derivative_1.Forward.Order2.coeffs = {-3./2., 2., -1./2.};

    Derivative_1.Forward.Order3.idxs = {0, 1, 2, 3};
    Derivative_1.Forward.Order3.coeffs = {-11./6., 3., -3./2., 1./3.};

    Derivative_1.Forward.Order4.idxs = {0, 1, 2, 3, 4};
    Derivative_1.Forward.Order4.coeffs = {-25./12., 4., -3., 4./3., -1./4.};

    Derivative_1.Forward.Order5.idxs = {0, 1, 2, 3, 4, 5};
    Derivative_1.Forward.Order5.coeffs = {-137./60., 5., -5., 10./3., -5./4., 1./5.};

    Derivative_1.Forward.Order6.idxs = {0, 1, 2, 3, 4, 5, 6};
    Derivative_1.Forward.Order6.coeffs = {-49./20., 6., -15./2., 20./3., -15./4., 6./5., -1./6.};

    // Central
    Derivative_1.Central.Order2.idxs = {-1, 1};
    Derivative_1.Central.Order2.coeffs = {-1./2., 1./2.};

    Derivative_1.Central.Order4.idxs = {-2, -1, 1, 2};
    Derivative_1.Central.Order4.coeffs = {1./12., -2./3., 2./3., -1./12.};

    Derivative_1.Central.Order6.idxs = {-3, -2, -1, 1, 2, 3};
    Derivative_1.Central.Order6.coeffs = {-1./60., 3./20., -3./4., 3./4., -3./20., 1./60.};

    Derivative_1.Central.Order8.idxs = {-4, -3, -2, -1, 1, 2, 3, 4};
    Derivative_1.Central.Order8.coeffs = {1./280., -4./105., 1./5., -4./5., 4./5., -1./5., 4./105., -1./280.};

    // Backward
    Derivative_1.Backward = Reverse (Derivative_1.Forward);


    // Derivée seconde
    // Forward
    Derivative_2.Forward.Order1.idxs = {0, 1, 2};
    Derivative_2.Forward.Order1.coeffs = {1., 2., -1./2.};

    Derivative_2.Forward.Order2.idxs = {0, 1, 2, 3};
    Derivative_2.Forward.Order2.coeffs = {2., -5., 4., -1.};

    Derivative_2.Forward.Order3.idxs = {0, 1, 2, 3, 4};
    Derivative_2.Forward.Order3.coeffs = {35./12., -26./3., 19./2., -14./3., 11./12.};

    Derivative_2.Forward.Order4.idxs = {0, 1, 2, 3, 4, 5};
    Derivative_2.Forward.Order4.coeffs = {15./4., -77./6., 107./6., -13., 61./12., -5./6.};

    Derivative_2.Forward.Order5.idxs = {0, 1, 2, 3, 4, 5, 6};
    Derivative_2.Forward.Order5.coeffs = {203./45., -87./5., 117./4., -254./9., 33./2., -27./5., 137./180.};

    Derivative_2.Forward.Order6.idxs = {0, 1, 2, 3, 4, 5, 6, 7};
    Derivative_2.Forward.Order6.coeffs = {469./90., -223./10., 879./20., -949./18., 41., -201./10., 1019./180., -7./10.};

    // Central
    Derivative_2.Central.Order2.idxs = {-1, 0, 1};
    Derivative_2.Central.Order2.coeffs = {1., -2.,  1.};

    Derivative_2.Central.Order4.idxs = {-2, -1, 0, 1, 2};
    Derivative_2.Central.Order4.coeffs = {-1./12., 4./3., -5./2., 4./3., -1./12.};

    Derivative_2.Central.Order6.idxs = {-3, -2, -1, 0, 1, 2, 3};
    Derivative_2.Central.Order6.coeffs = {1./90., -3./20., 3./2., -49./18., 3./2., -3./20., 1./90.};

    Derivative_2.Central.Order8.idxs = {-4, -3, -2, -1, 0, 1, 2, 3, 4};
    Derivative_2.Central.Order8.coeffs = {-1./560., 8./315., -1./5., 8./5., -205./72., 8./5., -1./5., 8./315., -1./560.};

    // Backward
    Derivative_2.Backward = Reverse (Derivative_2.Forward);
}


DFOrders Reverse(DFOrders O)
{
    DFOrders N;

    N.Order1 = Reverse (O.Order1);
    N.Order2 = Reverse (O.Order2);
    N.Order3 = Reverse (O.Order3);
    N.Order4 = Reverse (O.Order4);
    N.Order5 = Reverse (O.Order5);
    N.Order6 = Reverse (O.Order6);
    N.Order7 = Reverse (O.Order7);
    N.Order8 = Reverse (O.Order8);

    return N;
}

DFIdxCoeff Reverse (DFIdxCoeff D)
{
    DFIdxCoeff N = D;

    N.idxs = -1 * N.idxs;
    N.coeffs = -1. * N.coeffs;

    std::reverse(N.idxs.begin (), N.idxs.end ());
    std::reverse(N.coeffs.begin (), N.coeffs.end ());

    return N;
}

void DFOrderBuild (int der, ORDERS order, std::vector<int>* idxs, std::vector<double>* coeffs)
{
    DFStruct DF;

    if (der == 1)
    {
        switch (order)
        {
        case ORDER_1_FORWARD:
            *idxs = DF.Derivative_1.Forward.Order1.idxs;
            *coeffs = DF.Derivative_1.Forward.Order1.coeffs;
            break;
        case ORDER_1_CENTRAL:
            *idxs = DF.Derivative_1.Central.Order1.idxs;
            *coeffs = DF.Derivative_1.Central.Order1.coeffs;
            break;
        case ORDER_1_BACKWARD:
            *idxs = DF.Derivative_1.Backward.Order1.idxs;
            *coeffs = DF.Derivative_1.Backward.Order1.coeffs;
            break;
        case ORDER_2_FORWARD:
            *idxs = DF.Derivative_1.Forward.Order2.idxs;
            *coeffs = DF.Derivative_1.Forward.Order2.coeffs;
            break;
        case ORDER_2_CENTRAL:
            *idxs = DF.Derivative_1.Central.Order2.idxs;
            *coeffs = DF.Derivative_1.Central.Order2.coeffs;
            break;
        case ORDER_2_BACKWARD:
            *idxs = DF.Derivative_1.Backward.Order2.idxs;
            *coeffs = DF.Derivative_1.Backward.Order2.coeffs;
            break;
        case ORDER_3_FORWARD:
            *idxs = DF.Derivative_1.Forward.Order3.idxs;
            *coeffs = DF.Derivative_1.Forward.Order3.coeffs;
            break;
        case ORDER_3_CENTRAL:
            *idxs = DF.Derivative_1.Central.Order3.idxs;
            *coeffs = DF.Derivative_1.Central.Order3.coeffs;
            break;
        case ORDER_3_BACKWARD:
            *idxs = DF.Derivative_1.Backward.Order3.idxs;
            *coeffs = DF.Derivative_1.Backward.Order3.coeffs;
            break;
        case ORDER_4_FORWARD:
            *idxs = DF.Derivative_1.Forward.Order4.idxs;
            *coeffs = DF.Derivative_1.Forward.Order4.coeffs;
            break;
        case ORDER_4_CENTRAL:
            *idxs = DF.Derivative_1.Central.Order4.idxs;
            *coeffs = DF.Derivative_1.Central.Order4.coeffs;
            break;
        case ORDER_4_BACKWARD:
            *idxs = DF.Derivative_1.Backward.Order4.idxs;
            *coeffs = DF.Derivative_1.Backward.Order4.coeffs;
            break;
        case ORDER_5_FORWARD:
            *idxs = DF.Derivative_1.Forward.Order5.idxs;
            *coeffs = DF.Derivative_1.Forward.Order5.coeffs;
            break;
        case ORDER_5_CENTRAL:
            *idxs = DF.Derivative_1.Central.Order5.idxs;
            *coeffs = DF.Derivative_1.Central.Order5.coeffs;
            break;
        case ORDER_5_BACKWARD:
            *idxs = DF.Derivative_1.Backward.Order5.idxs;
            *coeffs = DF.Derivative_1.Backward.Order5.coeffs;
            break;
        case ORDER_6_FORWARD:
            *idxs = DF.Derivative_1.Forward.Order6.idxs;
            *coeffs = DF.Derivative_1.Forward.Order6.coeffs;
            break;
        case ORDER_6_CENTRAL:
            *idxs = DF.Derivative_1.Central.Order6.idxs;
            *coeffs = DF.Derivative_1.Central.Order6.coeffs;
            break;
        case ORDER_6_BACKWARD:
            *idxs = DF.Derivative_1.Backward.Order6.idxs;
            *coeffs = DF.Derivative_1.Backward.Order6.coeffs;
            break;
        case ORDER_7_FORWARD:
            *idxs = DF.Derivative_1.Forward.Order7.idxs;
            *coeffs = DF.Derivative_1.Forward.Order7.coeffs;
            break;
        case ORDER_7_CENTRAL:
            *idxs = DF.Derivative_1.Central.Order7.idxs;
            *coeffs = DF.Derivative_1.Central.Order7.coeffs;
            break;
        case ORDER_7_BACKWARD:
            *idxs = DF.Derivative_1.Backward.Order7.idxs;
            *coeffs = DF.Derivative_1.Backward.Order7.coeffs;
            break;
        case ORDER_8_FORWARD:
            *idxs = DF.Derivative_1.Forward.Order8.idxs;
            *coeffs = DF.Derivative_1.Forward.Order8.coeffs;
            break;
        case ORDER_8_CENTRAL:
            *idxs = DF.Derivative_1.Central.Order8.idxs;
            *coeffs = DF.Derivative_1.Central.Order8.coeffs;
            break;
        case ORDER_8_BACKWARD:
            *idxs = DF.Derivative_1.Backward.Order8.idxs;
            *coeffs = DF.Derivative_1.Backward.Order8.coeffs;
            break;
        }

    } else if (der == 2)
    {
        switch (order)
        {
        case ORDER_1_FORWARD:
            *idxs = DF.Derivative_2.Forward.Order1.idxs;
            *coeffs = DF.Derivative_2.Forward.Order1.coeffs;
            break;
        case ORDER_1_CENTRAL:
            *idxs = DF.Derivative_2.Central.Order1.idxs;
            *coeffs = DF.Derivative_2.Central.Order1.coeffs;
            break;
        case ORDER_1_BACKWARD:
            *idxs = DF.Derivative_2.Backward.Order1.idxs;
            *coeffs = DF.Derivative_2.Backward.Order1.coeffs;
            break;
        case ORDER_2_FORWARD:
            *idxs = DF.Derivative_2.Forward.Order2.idxs;
            *coeffs = DF.Derivative_2.Forward.Order2.coeffs;
            break;
        case ORDER_2_CENTRAL:
            *idxs = DF.Derivative_2.Central.Order2.idxs;
            *coeffs = DF.Derivative_2.Central.Order2.coeffs;
            break;
        case ORDER_2_BACKWARD:
            *idxs = DF.Derivative_2.Backward.Order2.idxs;
            *coeffs = DF.Derivative_2.Backward.Order2.coeffs;
            break;
        case ORDER_3_FORWARD:
            *idxs = DF.Derivative_2.Forward.Order3.idxs;
            *coeffs = DF.Derivative_2.Forward.Order3.coeffs;
            break;
        case ORDER_3_CENTRAL:
            *idxs = DF.Derivative_2.Central.Order3.idxs;
            *coeffs = DF.Derivative_2.Central.Order3.coeffs;
            break;
        case ORDER_3_BACKWARD:
            *idxs = DF.Derivative_2.Backward.Order3.idxs;
            *coeffs = DF.Derivative_2.Backward.Order3.coeffs;
            break;
        case ORDER_4_FORWARD:
            *idxs = DF.Derivative_2.Forward.Order4.idxs;
            *coeffs = DF.Derivative_2.Forward.Order4.coeffs;
            break;
        case ORDER_4_CENTRAL:
            *idxs = DF.Derivative_2.Central.Order4.idxs;
            *coeffs = DF.Derivative_2.Central.Order4.coeffs;
            break;
        case ORDER_4_BACKWARD:
            *idxs = DF.Derivative_2.Backward.Order4.idxs;
            *coeffs = DF.Derivative_2.Backward.Order4.coeffs;
            break;
        case ORDER_5_FORWARD:
            *idxs = DF.Derivative_2.Forward.Order5.idxs;
            *coeffs = DF.Derivative_2.Forward.Order5.coeffs;
            break;
        case ORDER_5_CENTRAL:
            *idxs = DF.Derivative_2.Central.Order5.idxs;
            *coeffs = DF.Derivative_2.Central.Order5.coeffs;
            break;
        case ORDER_5_BACKWARD:
            *idxs = DF.Derivative_2.Backward.Order5.idxs;
            *coeffs = DF.Derivative_2.Backward.Order5.coeffs;
            break;
        case ORDER_6_FORWARD:
            *idxs = DF.Derivative_2.Forward.Order6.idxs;
            *coeffs = DF.Derivative_2.Forward.Order6.coeffs;
            break;
        case ORDER_6_CENTRAL:
            *idxs = DF.Derivative_2.Central.Order6.idxs;
            *coeffs = DF.Derivative_2.Central.Order6.coeffs;
            break;
        case ORDER_6_BACKWARD:
            *idxs = DF.Derivative_2.Backward.Order6.idxs;
            *coeffs = DF.Derivative_2.Backward.Order6.coeffs;
            break;
        case ORDER_7_FORWARD:
            *idxs = DF.Derivative_2.Forward.Order7.idxs;
            *coeffs = DF.Derivative_2.Forward.Order7.coeffs;
            break;
        case ORDER_7_CENTRAL:
            *idxs = DF.Derivative_2.Central.Order7.idxs;
            *coeffs = DF.Derivative_2.Central.Order7.coeffs;
            break;
        case ORDER_7_BACKWARD:
            *idxs = DF.Derivative_2.Backward.Order7.idxs;
            *coeffs = DF.Derivative_2.Backward.Order7.coeffs;
            break;
        case ORDER_8_FORWARD:
            *idxs = DF.Derivative_2.Forward.Order8.idxs;
            *coeffs = DF.Derivative_2.Forward.Order8.coeffs;
            break;
        case ORDER_8_CENTRAL:
            *idxs = DF.Derivative_2.Central.Order8.idxs;
            *coeffs = DF.Derivative_2.Central.Order8.coeffs;
            break;
        case ORDER_8_BACKWARD:
            *idxs = DF.Derivative_2.Backward.Order8.idxs;
            *coeffs = DF.Derivative_2.Backward.Order8.coeffs;
            break;
        }

    }

    return;
}
