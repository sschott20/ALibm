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

#define SAMPLE_START 1.01
#define SAMPLE_END 1.04

using namespace std;
using namespace soplex;

double RangeReduction(float x)
{
    int exp;
    double sig = frexp(x, &exp);
    return sig;
}

double OutputCompensation(float x, double yp)
{
    int exp;
    double sig = frexp(x, &exp);
    return yp + exp;
}

double InverseOutputCompensation(float x, double yp)
{
    int exp;
    double sig = frexp(x, &exp);
    return yp - exp;
}

