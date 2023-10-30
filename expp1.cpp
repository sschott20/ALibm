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
#include "helper.hpp"

using namespace std;
using namespace soplex;

#define ROUND_2_INT(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))
#define HEX_TO_DOUBLE(h) ((doubleX){i : h}).d

typedef union
{
    double d;
    long long i;
} doubleX;

struct ReducedFloat
{
    int m;
    int j;
    double r1;
    double r2;
};

double Inv_L = HEX_TO_DOUBLE(0x40471547652B82FE);
double L1 = HEX_TO_DOUBLE(0x3F962E42FEF00000);
double L2 = HEX_TO_DOUBLE(0x3D8473DE6AF278ED);

// reduce x to the range [-ln2 / 64, ln2 / 64], X = (32M + J)ln2 / 64 + (R1 + R2)
// M, J : int
// R1, R2 : double
ReducedFloat
RangeReduction(float x)
{
    ReducedFloat rf;
    double r1, r2;
    int n, n1, n2;

    n = ROUND_2_INT(x * Inv_L);
    n2 = n % 32;

    printf("%d\n", n2);
    printf("%d\n", n);

    n1 = n - n2;

    if (abs(n) > 2 ^ 9)
    {
        r1 = (x - n1 * L1) - n2 * L1;
    }
    else
    {
        r1 = x - n * (double)L1;
    }

    r2 = -n * (double)L2;

    rf.j = n2;
    rf.m = n1 / 32;

    rf.r1 = r1;
    rf.r2 = r2;

    cout << rf.j << " " << rf.m << " " << rf.r1 << " " << rf.r2 << endl;
    cout << (double)(32 * rf.m + rf.j) * ((double)L1 + (double)L2) + (double)(rf.r1 + rf.r2) << endl;
    return rf;
}

int main(int argc, char *argv[])
{
    ReducedFloat rf = RangeReduction(1.125);

    return 0;
}