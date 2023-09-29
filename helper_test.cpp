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

using namespace std;
using namespace soplex;
using half_float::half;

struct RndInterval
{
    half x_orig;
    half y;
    double l;
    double u;

    double x_rr;
    double yp;
    double lp;
    double up;
};

struct Polynomial
{
    int termsize;
    vector<double> coefficients;
};

void RangeReduction(half x)
{
    int exp;
    double sig = half_float::frexp(x, &exp);
    // exp -= 1
    // sig * 2;
    // print_binary_half(x);
    // std::cout << x << " makes " << sig << " * 2^" << exp - 1 << '\n';

    // return sig * 2;
    printf("x = %4.15f sig = %4.15f exp = %d\n", (double)x, sig, exp);
}

double OutputCompensation(half x, double yp)
{
    int exp;
    double sig = half_float::frexp(x, &exp);
    // sig *= 2

    // return yp + exp - 1;
    // printf("x = %4.15f yp = %4.15f, exp = %d\n", (double)x, yp, exp);

    return yp + exp;
}

double InverseOutputCompensation(half x, double yp)
{
    int exp;
    double sig = half_float::frexp(x, &exp);

    // return yp - (exp + 1);
    // printf("x = %4.15f yp = %4.15f, exp = %d ret = %4.15f\n", (double)x, yp, exp, yp - exp + 1);
    return yp - exp;
}

int main(int argc, char *argv[])
{
    if (argc == 2)
    {
        half input((double)atof(argv[1]));
        RangeReduction(input);
    }
    else if (argc == 3)
    {
        double in = atof(argv[1]);
        double yp = atof(argv[2]);
        half x(in);
        printf("OC = %4.30f Oc-1 = %4.30f\n", OutputCompensation(x, yp), InverseOutputCompensation(x, yp));

        ;
    }

    return 0;
}
