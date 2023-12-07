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

void print_binary_half(half h);
vector<RndInterval> GenerateFloatSample(int sample_size);
double EvaulutePoly(Polynomial P, double xval);
vector<RndInterval> Verify(vector<RndInterval> L2, Polynomial P);
vector<RndInterval> CalcRndIntervals(vector<RndInterval> X);
vector<RndInterval> CalcRedIntervals(vector<RndInterval> X);
Polynomial GeneratePolynomial(vector<RndInterval> L);
