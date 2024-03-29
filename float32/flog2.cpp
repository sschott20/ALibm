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
#include <math.h>
#include <bitset>
#include <iomanip>
#include "luts.h"
#include "fhelper.hpp"

using namespace std;
using namespace soplex;

int ComputeSpecialCase(float x)
{
    floatX fx;
    fx.f = x;
    if (x <= 0.0)
    {
        return -1;
    }
    else if (fx.x >= 0x7F800000)
    {
        return -1;
    }
    else if (fx.f == 2.0037834644317626953125000000000000000000)
    {
        return -1;
    }

    return 0;
}
double RangeReduction(float x)
{
    // reduces x to [0.5, 1)
    int exp;
    double sig = (double)frexp(x, &exp);
    return sig;
    // floatX fix, fit;
    // int m = 0;
    // fix.f = x;
    // if (fix.x < 0x800000)
    // {
    //     fix.f *= pow(2, 23);
    //     m -= 23;
    // }
    // m += fix.x >> 23;
    // m -= 127;
    // fix.x &= 0x007FFFFF;
    // fix.x |= 0x3F800000;

    // fit.x = fix.x & 0x007F0000;
    // int FIndex = fit.x >> 16;
    // fit.x |= 0x3F800000;
    // double F = fit.f;

    // double f = fix.f - F;
    // return f * log2OneByF[FIndex];
}

double OutputCompensation(float x, double yp)
{
    int exp;
    double sig = (double)frexp(x, &exp);
    return yp + exp;
    // floatX fix, fit;

    // int m = 0;
    // fix.f = x;
    // if (fix.x < 0x800000)
    // {
    //     fix.f *= pow(2, 23);
    //     m -= 23;
    // }
    // m += fix.x >> 23;
    // m -= 127;
    // fix.x &= 0x007FFFFF;
    // fix.x |= 0x3F800000;

    // fit.x = fix.x & 0x007F0000;
    // int FIndex = fit.x >> 16;

    // return yp + log2Lut[FIndex] + m;
}

double GuessInitial(float x, double &lb, double &ub)
{
    // int exp;
    // double sig = (double)frexp(x, &exp);
    // return yp - exp;
    // lb = log1p(xp) / log(2);
    // ub = log1p(xp) / log(2);
    double xrr = RangeReduction(x);
    double yp = log2(xrr);
    lb = yp;
    ub = yp;
    return 0;
}

float EvaluateFunction(mpfr_t y, double x)
{
    mpfr_set_d(y, x, MPFR_RNDN);
    mpfr_log2(y, y, MPFR_RNDN);
    float h = mpfr_get_flt(y, MPFR_RNDN);
    return h;
}

#define GROW 10
int main()
{

    printf("Generating FloatSample...\n");
    vector<RndInterval> X = GenerateFloatSample(10, 1, 2);

    printf("Generating all float values...\n");

    vector<RndInterval> Test = GenerateFloatSample(-1, 1, 2);
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

        Incorrect = Verify(L2, P, 1);
        if (Incorrect.size() > 0)
        {
            printf("FAILED IN SAMPLE: %ld\n", Incorrect.size());

            // return 0;
        }
        Incorrect = Verify(Test, P, 0);
        print_poly(P);
        printf("----------------------------------------\n");

    } while (Incorrect.size() > 0);
    // Test = GenerateFloatSample(-1, 2, 10);

    // Test = GenerateFloatSample(-1, 0.1, 99999999999);
    FullTest(P, 2, 20);
    print_poly(P);
    printf("Finished!\n");
}