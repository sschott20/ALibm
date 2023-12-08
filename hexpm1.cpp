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
#include "luts.h"
#include "hhelper.hpp"

using namespace std;
using namespace soplex;
using half_float::half;

int ComputeSpecialCase(half x)
{
    mpfr_t y_mpfr;
    mpfr_init2(y_mpfr, 200);
    half y = EvaluateFunction(y_mpfr, x);
    // check if y is infinity
    if (y == INFINITY)
    {
        return -1;
    }
    if (y == -1)
    {
        return -1;
    }
    return 0;
}

double RangeReduction(half x)
{
    // double xp = x * 92.332482616893656768297660164535045623779296875;
    // int N = (int)xp;
    // return x - N * 0.01083042469624914509729318723429969395510852336883544921875;
    return (double)x;
}

double OutputCompensation(half x, double yp)
{
    // double xp = x * 92.332482616893656768297660164535045623779296875;
    // int N = (int)xp;
    // int N2 = N % 64;
    // if (N2 < 0)
    //     N2 += 64;
    // int N1 = N - N2;
    // int M = N1 / 64;

    // return yp * ldexp(exp2JBy64[N2], M);
    return yp;
}

double InverseOutputCompensation(half x, double yp)
{
    // double xp = x * 92.332482616893656768297660164535045623779296875;
    // int N = (int)xp;
    // int N2 = N % 64;
    // if (N2 < 0)
    //     N2 += 64;
    // int N1 = N - N2;
    // int M = N1 / 64;
    // int J = N2;
    // return (2 << M) * ((2 << (J / 64)) + (2 << J / 64) * yp) - 1;
    return yp;
}

half EvaluateFunction(mpfr_t y, double x)
{
    mpfr_set_d(y, x, MPFR_RNDN);
    mpfr_expm1(y, y, MPFR_RNDN);
    half h = (half)mpfr_get_d(y, MPFR_RNDN);
    return h;
}

int main()
{
    printf("Generating FloatSample...\n");
    vector<RndInterval> X = GenerateFloatSample(10, -10, 10);

    printf("Generating all float values...\n");
    vector<RndInterval> Test = GenerateFloatSample(-1, -1, 1);
    vector<RndInterval> Incorrect;
    Polynomial P;

    do
    {
        if (Incorrect.size() < 10)
        {
            for (int i = 0; i < Incorrect.size(); i++)
            {
                RndInterval I = Incorrect.at(i);
                X.push_back(I);
            }
        }
        else
        {
            for (int i = 0; i < 10; i++)
            {
                RndInterval I = Incorrect.at(floor(i * Incorrect.size() / 11));
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

    for (int i = 0; i < P.termsize; i++)
    {
        printf("%5.60f\n", P.coefficients.at(i));
    }

    FILE *fptr = fopen("dump/poly.txt", "w");
    for (int i = 0; i < P.termsize; i++)
    {
        fprintf(fptr, "%5.60f\n", P.coefficients.at(i));
    }
    fclose(fptr);

    printf("Finished!\n");
}