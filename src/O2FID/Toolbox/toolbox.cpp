#include "toolbox.h"

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
