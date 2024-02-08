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

double RangeReduction(half x);
double OutputCompensation(half x, double yp);
double InverseOutputCompensation(half x, double yp);
half EvaluateFunction(mpfr_t y, double x);

void print_binary_half(half h);

vector<RndInterval> GenerateFloatSample(int sample_size, float min, float max);

vector<RndInterval> CalcRndIntervals(vector<RndInterval> X);

vector<RndInterval> CalcRedIntervals(vector<RndInterval> X);

Polynomial GeneratePolynomial(vector<RndInterval> L);

double EvaulutePoly(Polynomial P, double xval);

vector<RndInterval> Verify(vector<RndInterval> L2, Polynomial P);

int ComputeSpecialCase(half x);

void print_poly(Polynomial P);
void print_binary_half(half h);
