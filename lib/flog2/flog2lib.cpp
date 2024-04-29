#include "flog2lib.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <mpfr.h>
#include <cassert>
using namespace std;

Polynomial PTest = InitPoly(0);

std::vector<Polynomial> flog2_poly = {
    InitPoly(0),
    InitPoly(1),
    InitPoly(2),
    InitPoly(3),
    InitPoly(4)

};

Polynomial InitPoly(int index)

{
    Polynomial P;
    P.termsize = 0;
    for (int i = 0; i < 7; i++)
    {
        doubleX dx;
        dx.x = flog2_coef[index][i];
        P.coefficients.push_back(dx.d);
    }
    return P;
}
double EvaulutePoly(int index, double xval)
{
    double acc = 0.0;
    double power = 1.0;
    doubleX dx;
    for (int i = 0; i < 7; i++)
    {
        dx.x = flog2_coef[index][i];
        acc += dx.d * power;
        power *= xval;
    }

    return acc;
}
double RangeReduction(float x)
{
    // return x;

    int exp;
    double sig = (double)frexp(x, &exp);
    return sig;
}

double OutputCompensation(float x, double yp)
{
    int exp;
    double sig = (double)frexp(x, &exp);
    return yp + exp;
}
float EvaluateFunction(mpfr_t y, float x)
{
    mpfr_set_d(y, x, MPFR_RNDN);
    mpfr_log2(y, y, MPFR_RNDN);
    float h = mpfr_get_flt(y, MPFR_RNDN);
    return h;
}

int SelectPoly(float x)
{
    if (x < 2.0)
    {
        return -1;
    }
    else if (x < 2.2)
    {
        return 0;
    }
    else if (x < 2.4)
    {
        return 1;
    }
    else if (x < 2.6)
    {
        return 2;
    }
    else if (x < 2.8)
    {
        return 3;
    }
    else if (x < 3.0)
    {
        return 4;
    }
    else
    {
        return -1;
    }
}
float flog2(float x)
{
    int poly_ind = SelectPoly(x);
    double range_reduced = RangeReduction(x);
    //  double eval = EvaulutePoly(flog2_poly[poly_ind], range_reduced);
    // Polynomial P = InitPoly(poly_ind);

    double eval = EvaulutePoly(poly_ind, range_reduced);
    double oc = OutputCompensation(x, eval);
    float y = (float)oc;
    return y;
}

#define LOW 2.0
#define HIGH 3.0

int main()
{
    // for (int i = 0; i < flog2_poly.size(); i++)
    // {
    //     for (int j = 0; j < flog2_poly[i].coefficients.size(); j++)
    //     {
    //         cout << flog2_poly[i].coefficients[j] << endl;
    //     }
    // }
    size_t n = 0;
    size_t incorrect = 0;
    for (float x = LOW; x < HIGH; x = nextafterf(x, numeric_limits<float>::infinity()))
    {
        n++;
        float y = flog2(x);

        mpfr_t y_mpfr;
        mpfr_init2(y_mpfr, 500);
        float y2 = EvaluateFunction(y_mpfr, x);
        if (y != y2)
        {
            incorrect++;
            cout << "Mismatch detected for x = " << x << endl;
            cout << "Expected: " << y << endl;
            cout << "Actual: " << y2 << endl;

            exit(0);
        }
    }
    cout << "Total number of tests: " << n << endl;
    cout << "Number of incorrect results: " << incorrect << endl;
    cout << "Test complete." << endl;
    return 0;
}