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

typedef union
{
    double d;
    long long i;
} doubleX;

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
vector<RndInterval> GenerateFloatSample(int sample_size)
{
    vector<RndInterval> X;
    size_t n = 0;
    half h;
    half zero(0.0);
    half inf((0x7BFF));
    if (sample_size == -1)
    {
        for (h = half_float::nextafter(zero, inf); h != inf; h = half_float::nextafter(h, inf))
        {
            n++;
            RndInterval I;
            I.x_orig = h;
            I.x_rr = RangeReduction(h);
            X.push_back(I);
        }
    }
    else
    {
        int skip = (1 << 16) / sample_size;
        for (h = half_float::nextafter(zero, inf); h != inf; h = half_float::nextafter(h, inf))
        {
            n++;
            RndInterval I;
            I.x_orig = h;
            I.x_rr = RangeReduction(h);
            X.push_back(I);
            for (int i = 0; i < skip; i++)
            {
                h = half_float::nextafter(h, inf);
            }
        }
    }
    // printf("sample size = %ld\n", n);
    return X;
}

vector<RndInterval> CalcRndIntervals(vector<RndInterval> X)
{
    FILE *fptr;
    fptr = fopen("dump/RndIntervalHalf.txt", "w");

    half y;
    double l, u;

    mpfr_t y_mpfr, x_mpfr;
    mpfr_inits2(200, y_mpfr, x_mpfr, NULL);

    vector<RndInterval> L;

    for (size_t i = 0; i < X.size(); i++)
    {

        RndInterval I;

        I.x_orig = X.at(i).x_orig;
        I.x_rr = X.at(i).x_rr;
        y = EvaluateFunction(y_mpfr, (double)X.at(i).x_orig);

        half inf(half_float::half(0x7BFF));

        half lhalf(half_float::nextafter((half)y, -inf));
        half uhalf(half_float::nextafter((half)y, inf));

        l = (y + (double)lhalf) / 2;
        u = (y + (double)uhalf) / 2;

        while ((half)l != (half)y)
        {

            l = nextafterf(l, INFINITY);
        }
        assert((half)l == (half)y);

        while ((half)u != (half)y)
        {

            u = nextafterf(u, -INFINITY);
        }
        assert((half)u == (half)y);
        assert(l < u);

        fprintf(fptr, "%4.15f %4.15f %4.40f %4.40f \n", (double)X.at(i).x_orig, (double)y, l, u);
        I.y = y;
        I.l = l;
        I.u = u;

        L.push_back(I);
    }
    return L;
}

vector<RndInterval> CalcRedIntervals(vector<RndInterval> X)
{
    FILE *fptr;
    fptr = fopen("dump/RedIntervalHalf.txt", "w");
    half yp;
    double lp, up;

    mpfr_t yp_mpfr, xrr_mpfr;
    mpfr_inits2(200, yp_mpfr, xrr_mpfr, NULL);

    vector<RndInterval> L;

    for (size_t i = 0; i < X.size(); i++)
    {
        RndInterval I;

        I.x_orig = X.at(i).x_orig;
        I.x_rr = X.at(i).x_rr;

        I.y = X.at(i).y;

        I.l = X.at(i).l;
        I.u = X.at(i).u;

        lp = InverseOutputCompensation(X.at(i).x_orig, X.at(i).l);
        up = InverseOutputCompensation(X.at(i).x_orig, X.at(i).u);

        while (OutputCompensation(X.at(i).x_orig, lp) < X.at(i).l)
        {
            lp = nextafterf(lp, INFINITY);
        }

        while (OutputCompensation(X.at(i).x_orig, up) > X.at(i).u)
        {
            up = nextafterf(up, -INFINITY);
        }

        assert(lp <= up);
        fprintf(fptr, "%4.15f %4.40f %4.40f %4.40f \n", (double)X.at(i).x_orig, X.at(i).x_rr, lp, up);
        I.lp = lp;
        I.up = up;

        L.push_back(I);
    }
    return L;
}

Polynomial GeneratePolynomial(vector<RndInterval> L)
{

    for (int termsize = 1; termsize < 30; termsize++)
    {
        SoPlex mysoplex;
        mysoplex.setBoolParam(SoPlex::RATFACJUMP, true);
        mysoplex.setIntParam(SoPlex::SOLVEMODE, 2);
        mysoplex.setIntParam(SoPlex::CHECKMODE, 2);
        mysoplex.setIntParam(SoPlex::SYNCMODE, 1);
        mysoplex.setIntParam(SoPlex::READMODE, 1);
        mysoplex.setRealParam(SoPlex::FEASTOL, 0.0);
        mysoplex.setRealParam(SoPlex::OPTTOL, 0.0);
        mysoplex.setRealParam(SoPlex::EPSILON_ZERO, 0.0);
        mysoplex.setRealParam(SoPlex::EPSILON_FACTORIZATION, 0.0);
        mysoplex.setRealParam(SoPlex::EPSILON_UPDATE, 0.0);
        mysoplex.setRealParam(SoPlex::EPSILON_PIVOT, 0.0);
        mysoplex.setIntParam(SoPlex::VERBOSITY, 0);
        mysoplex.setRealParam(SoPlex::TIMELIMIT, 5 * 60);
        mysoplex.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
        DSVectorRational dummycol(0);

        for (int i = 0; i < termsize; i++)
        {
            auto column = LPColRational(1.0, dummycol, infinity, -infinity);
            mysoplex.addColRational(column);
        }

        for (int i = 0; i < L.size(); i++)
        {
            DSVectorRational row1(termsize);
            Rational acc(1.0);
            Rational xval(L.at(i).x_rr);
            row1.add(0, 1.0);

            for (int j = 1; j < termsize; j++)
            {
                acc = acc * xval;
                row1.add(j, acc);
            }
            double lbnd = (L.at(i).lp);
            double ubnd = (L.at(i).up);

            mysoplex.addRowRational(LPRowRational(lbnd, row1, ubnd));
        }

        mysoplex.writeFileRational("dump/dump_rational.lp", NULL, NULL, NULL);
        mysoplex.writeFileReal("dump/dump_real.lp", NULL, NULL, NULL);

        SPxSolver::Status stat;
        cout << "Solving... " << termsize << "\n";

        stat = mysoplex.optimize();
        Polynomial P;

        if (stat == SPxSolver::OPTIMAL)
        {
            DVectorRational prim(termsize);
            mysoplex.getPrimalRational(prim);

            P.termsize = termsize;
            for (int i = 0; i < termsize; i++)
            {
                P.coefficients.push_back(prim[i]);
            }

            return P;
        }
        else if (stat == SPxSolver::UNBOUNDED)
        {
            DVectorRational prim(termsize);
            mysoplex.getPrimalRational(prim);

            P.termsize = termsize;
            std::cout << "Unbounded solution\n";
            for (int i = 0; i < termsize; i++)
            {
                P.coefficients.push_back(0.0);
            }
        }
        std::cout << "Status: " << stat << std::endl;
    }

    Polynomial N;
    return N;
}

double EvaulutePoly(Polynomial P, double xval)
{
    double acc = 0.0;
    double power = 1.0;

    for (int i = 0; i < P.termsize; i++)
    {

        acc += P.coefficients.at(i) * power;
        power *= xval;
    }

    return acc;
}

vector<RndInterval> Verify(vector<RndInterval> L2, Polynomial P)
{
    size_t n = 0;
    size_t correct = 0;
    half inf(half_float::half(0x7BFF));
    vector<RndInterval> Incorrect;
    for (size_t i = 0; i < L2.size(); i++)
    {
        n++;
        half h(L2.at(i).x_orig);

        double x = (double)h;
        double range_reduced = RangeReduction(h);
        double eval = EvaulutePoly(P, range_reduced);

        half y(OutputCompensation(h, eval));

        mpfr_t y_mpfr;
        mpfr_inits2(200, y_mpfr, NULL);
        half oracle = EvaluateFunction(y_mpfr, x);

        if (y != oracle)
        {
            RndInterval I;
            I.x_orig = h;
            I.x_rr = range_reduced;
            Incorrect.push_back(I);
        }
        else
        {
            correct += 1;
        }
    }

    printf("In sample, Correct: %ld Incorrect: %ld\n", correct, n - correct);
    return Incorrect;
}