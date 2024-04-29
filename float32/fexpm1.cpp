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
    // return (double)x;
    double xp = x * 92.332482616893656768297660164535045623779296875;
    int N = (int)xp;

    return x - N *
                   0.01083042469624914509729318723429969395510852336883544921875;
}

double OutputCompensation(float x, double yp)
{
    // return yp;
    double xp = x * 92.332482616893656768297660164535045623779296875;
    int N = (int)xp;
    int N2 = N % 64;
    if (N2 < 0)
        N2 += 64;
    int N1 = N - N2;
    int M = N1 / 64;

    return yp * ldexp(exp2JBy64[N2], M);
}

double GuessInitial(float x, double &lb, double &ub)
{

    double xrr = RangeReduction(x);
    double yp = expm1(xrr);
    lb = yp;
    ub = yp;
    return 0;
}

float EvaluateFunction(mpfr_t y, float x)
{
    mpfr_set_flt(y, x, MPFR_RNDN);
    mpfr_expm1(y, y, MPFR_RNDN);
    float h = mpfr_get_flt(y, MPFR_RNDN);
    return h;
}

#define GROW 20
#define SPACING 0.1
#define OVERLAP 0.00001
int main()
{
    for (float low = 2; low < 3; low += SPACING)
    {
        float high = low + SPACING + OVERLAP;
        // printf("Generating FloatSample...\n");
        vector<RndInterval> X = GenerateFloatSample(1, low, high);

        // printf("Generating all float values...\n");

        vector<RndInterval> Test = GenerateFloatSample(-1, low, high);
        vector<RndInterval> Incorrect;
        Polynomial P;

        do
        {
            printf("Low: %f, High: %f\n", low, high);
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
                    for (size_t j = 0; j < X.size(); j++)
                    {
                        if (X.at(j).x_orig == I.x_orig || X.at(j).x_rr == I.x_rr)
                        {
                            printf("Duplicate\n");
                            break;
                        }
                    }
                    X.push_back(I);
                }
            }
            printf("Sample size: %ld\n", X.size());
            // printf("Generating RndIntervals...\n");
            vector<RndInterval> L = CalcRndIntervals(X);

            // printf("Generating RedIntervals...\n");
            vector<RndInterval> L2 = CalcRedIntervals(L);

            // printf("Generating Polynomial...\n");
            P = GeneratePolynomial(L2);

            Incorrect = Verify(L2, P, 1);
            if (Incorrect.size() > 0)
            {
                printf("FAILED IN SAMPLE: %ld\n", Incorrect.size());
                // return -1;
            }
            Incorrect = Verify(Test, P, 0);
            // print_poly(P);
            printf("----------------------------------------\n");

        } while (Incorrect.size() > 0);

        FullTest(P, low, high);
        std::string filename = "fexpm1poly/fexmp1_" + std::to_string(low) + "_" + std::to_string(high) + ".txt";
        writePolynomialToFile(P, filename);
        printf("Successfully written to file: %s\n", filename.c_str());
        print_poly(P);
    }
    printf("Finished!\n");
}