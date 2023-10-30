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

typedef union
{
    float f;
    int i;
} floatX;

struct Polynomial
{
    int termsize;
    vector<double> coefficients;
};

int write_rnd_interval_to_file(vector<RndInterval> X, const char *filename)
{
    FILE *fptr;
    fptr = fopen(filename, "wb");

    for (size_t i = 0; i < X.size(); i++)
    {

        fwrite(&X.at(i), sizeof(RndInterval), 1, fptr);
    }
    fclose(fptr);

    return 0;
}
vector<RndInterval> read_rnd_interval_from_file(const char *filename)
{
    vector<RndInterval> X;
    FILE *fptr;
    fptr = fopen(filename, "rb");
    RndInterval x1;

    while (fread(&x1, sizeof(RndInterval), 1, fptr) == 1)
    {
        X.push_back(x1);
    }
    fclose(fptr);

    return X;
}