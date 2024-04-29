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

    bool operator<(RndInterval &other) const
    {
        // return std::make_tuple(fromX, fromY, toX, toY) < std::make_tuple(other.fromX, other.fromY, other.toX, other.toY);
        return x_orig < other.x_orig;
    }

    bool operator==(RndInterval &other) const
    {
        // declare how 2 variable of type posToMove should be compared with ==
        // return std::make_tuple(fromX, fromY, toX, toY) == std::make_tuple(other.fromX, other.fromY, other.toX, other.toY);
        return x_orig == other.x_orig;
    }
};

struct Polynomial
{
    int termsize;
    vector<double> coefficients;
};
typedef union
{
    float f;
    unsigned x;
} floatX;

typedef union
{
    double d;
    unsigned long long int x;
} doubleX;
void FullTest(Polynomial P, float min = 0.0, float max = INFINITY);
double RangeReduction(float x);
double OutputCompensation(float x, double yp);
double InverseOutputCompensation(float x, double yp);
float EvaluateFunction(mpfr_t y, float x);
vector<RndInterval> GenerateFloatSample(int sample_size, float min, float max);
vector<RndInterval> CalcRndIntervals(vector<RndInterval> X);
vector<RndInterval> CalcRedIntervals(vector<RndInterval> X);
Polynomial GeneratePolynomial(vector<RndInterval> L);
long double EvaulutePoly(Polynomial P, double xval);
vector<RndInterval> Verify(vector<RndInterval> L2, Polynomial P, int debug);
vector<RndInterval> VerifyAdd(Polynomial P, long last_size, float min, float max);
int ComputeSpecialCase(float x);
void print_poly(Polynomial P);
float generateRandomFloat(float a, float b);
double GuessInitial(float x, double &lb, double &ub);
bool CompareX(RndInterval a, RndInterval b);
void writePolynomialToFile(const Polynomial &poly, const std::string &filename);
Polynomial readPolynomialFromFile(const std::string &filename);
