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

std::vector<double> Order (std::vector<double> err, std::vector<double> h)
{
    size_t N = std::min (err.size (), h.size ());

    std::vector<double> order (N, 0.);

    for (size_t i = 1; i < N; ++i)
        order.at (i) = (std::log (err.at (i)) - std::log( err.at (i-1))) / (std::log (h.at (i)) - std::log(h.at (i-1)));

    return order;
}
