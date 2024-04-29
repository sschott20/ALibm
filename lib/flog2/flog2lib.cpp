#include "flog2lib.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <mpfr.h>
#include <cassert>
#include <time.h>
using namespace std;

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
    double eval = EvaulutePoly(poly_ind, range_reduced);
    double oc = OutputCompensation(x, eval);
    float y = (float)oc;
    return y;
}

#define LOW 2.0
#define HIGH 3.0

int main()
{

    size_t n = 0;
    size_t incorrect = 0;

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
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

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    n = 0;
    double delta_ms = (((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000) / 1000);
    cout << "ALibm" << endl;
    cout << "Total number of tests: " << n << endl;
    cout << "Total time taken: " << delta_ms << " ms" << endl;
    cout << "Number of incorrect results: " << incorrect << endl;
    cout << "Test complete." << endl;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (float x = LOW; x < HIGH; x = nextafterf(x, numeric_limits<float>::infinity()))
    {
        n++;
        float y = log2(x);
        mpfr_t y_mpfr;
        mpfr_init2(y_mpfr, 500);
        float y2 = EvaluateFunction(y_mpfr, x);
        if (y != y2)
        {
            incorrect++;
            // cout << "Mismatch detected for x = " << x << endl;
            // cout << "Expected: " << y << endl;
            // cout << "Actual: " << y2 << endl;
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    delta_ms = (((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000) / 1000);
    cout << "C Std lib" << endl;
    cout << "Total number of tests: " << n << endl;
    cout << "Total time taken: " << delta_ms << " ms" << endl;
    cout << "Number of incorrect results: " << incorrect << endl;
    cout << endl
         << endl;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    for (int i = 0; i < 10; i++)
    {
        for (float x = LOW; x < HIGH; x = nextafterf(x, numeric_limits<float>::infinity()))
        {
            float y = log2(x);
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    delta_ms = (((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000) / 1000);
    cout << "ALibm" << endl;
    cout << "Total time taken: " << delta_ms << " ms" << endl;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    for (int i = 0; i < 10; i++)
    {
        for (float x = LOW; x < HIGH; x = nextafterf(x, numeric_limits<float>::infinity()))
        {
            float y = flog2(x);
        }
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    delta_ms = (((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000) / 1000);
    cout << "C Std Lib" << endl;
    cout << "Total time taken: " << delta_ms << " ms" << endl;

    return 0;
}