#include <stdio.h>
#include <iostream>
#include <gmp.h>
#include <mpfr.h>
#include <inttypes.h>
#include <math.h>
#include <bits/stdc++.h>
#include <random>
#include <soplex.h>
#include <cstdlib>
#include "half.hpp"
#include <math.h>
#include <bitset>
#include <iomanip>
#include "hhelper.hpp"

using namespace std;
using namespace soplex;
using half_float::half;

int ComputeSpecialCase(half x)
{
    mpfr_t y_mpfr;
    mpfr_init2(y_mpfr, 200);
    half y = EvaluateFunction(y_mpfr, x);

    if (x <= 0)
    {
        return -1;
    }

    return 0;
}
double RangeReduction(half x)
{
    // reduces x to [0.5, 1)
    int exp;
    double sig = half_float::frexp(x, &exp);
    return sig;
}

double OutputCompensation(half x, double yp)
{
    int exp;
    double sig = half_float::frexp(x, &exp);
    return yp + exp;
}

double InverseOutputCompensation(half x, double yp)
{
    int exp;
    double sig = half_float::frexp(x, &exp);
    return yp - exp;
}

half EvaluateFunction(mpfr_t y, double x)
{
    mpfr_set_d(y, x, MPFR_RNDN);
    mpfr_log2(y, y, MPFR_RNDN);
    half h = (half)mpfr_get_d(y, MPFR_RNDN);
    return h;
}

#define GROW 10
int main()
{
    printf("Generating FloatSample...\n");
    vector<RndInterval> X = GenerateFloatSample(100, -INFINITY, INFINITY);

    printf("Generating all float values...\n");
    vector<RndInterval> Test = GenerateFloatSample(-1, 0, -1);
    vector<RndInterval> Incorrect;
    Polynomial P;

    do
    {
        if (Incorrect.size() < GROW)
        {
            for (int i = 0; i < Incorrect.size(); i++)
            {
                RndInterval I = Incorrect.at(i);
                X.push_back(I);
            }
        }
        else
        {
            for (int i = 0; i < GROW; i++)
            {
                RndInterval I = Incorrect.at(floor(i * Incorrect.size() / GROW));
                X.push_back(I);
            }
        }
        printf("Sample size: %ld\n", X.size());
        printf("Generating RndIntervals...\n");
        vector<RndInterval> L = CalcRndIntervals(X);

        printf("Generating RedIntervals...\n");
        vector<RndInterval> L2 = CalcRedIntervals(L);

        printf("Generating Polynomial...\n");
        P = GeneratePolynomial(L2);

        Incorrect = Verify(L2, P);
        if (Incorrect.size() > 0)
        {
            printf("FAILED IN SAMPLE: %ld\n", Incorrect.size());
            return 0;
        }
        Incorrect = Verify(Test, P);
    } while (Incorrect.size() > 0);

    print_poly(P);
    printf("Finished!\n");
}