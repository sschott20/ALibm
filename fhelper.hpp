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

using namespace std;
using namespace soplex;
struct RndInterval
{
    float x_orig;
    float y;
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

double RangeReduction(float x);
double OutputCompensation(float x, double yp);
double InverseOutputCompensation(float x, double yp);
float EvaluateFunction(mpfr_t y, double x);
vector<RndInterval> GenerateFloatSample(int sample_size, float min, float max);
vector<RndInterval> CalcRndIntervals(vector<RndInterval> X);
vector<RndInterval> CalcRedIntervals(vector<RndInterval> X);
Polynomial GeneratePolynomial(vector<RndInterval> L);
double EvaulutePoly(Polynomial P, double xval);
vector<RndInterval> Verify(vector<RndInterval> L2, Polynomial P, int debug);
int ComputeSpecialCase(float x);
void print_poly(Polynomial P);
float generateRandomFloat(float a, float b);
