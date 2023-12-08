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

void print_binary_half(half h)
{
    short n = *(short *)&h;

    unsigned i;
    for (i = 1 << 15; i > 0; i = i / 2)
        (n & i) ? printf("1") : printf("0");
    printf("\n");
}

double RangeReduction(half x)
{
    // reduces x to [0.5, 1)
    int exp;
    double sig = half_float::frexp(x, &exp);
    return sig;
}

double OutputCompensation(half x, double yp)
{
    int exp;
    double sig = half_float::frexp(x, &exp);
    return yp + exp;
}

double InverseOutputCompensation(half x, double yp)
{
    int exp;
    double sig = half_float::frexp(x, &exp);
    return yp - exp;
}

half EvaluateFunction(mpfr_t y, double x)
{
    mpfr_set_d(y, x, MPFR_RNDN);
    mpfr_log2(y, y, MPFR_RNDN);
    half h = (half)mpfr_get_d(y, MPFR_RNDN);
    return h;
}

vector<RndInterval> GenerateFloatSample(int sample_size, float min, float max)
{
    vector<RndInterval> X;
    size_t n = 0;
    half h;
    half start(min);
    half end(max);
    half inf((0x7BFF));

    if (max == -1)
    {
        end = inf;
    }

    if (max = INFINITY)
    {
        end = inf;
    }
    if (sample_size == -1)
    {
        for (h = half_float::nextafter(start, end); h < end; h = half_float::nextafter(h, end))
        {
            if (h == inf || h == -inf)
            {
                break;
            }
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

        for (h = half_float::nextafter(start, end); h < end; h = half_float::nextafter(h, end))
        {
            if (h == inf || h == -inf)
            {
                break;
            }
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

        printf("x: %4.15f y: %4.15f l: %4.15f u: %4.15f\n", (double)X.at(i).x_orig, (double)y, l, u);
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

int main()
{
    printf("Generating FloatSample...\n");
    vector<RndInterval> X = GenerateFloatSample(100, 0, -1);

    printf("Generating all float values...\n");
    vector<RndInterval> Test = GenerateFloatSample(-1, 0, -1);
    vector<RndInterval> Incorrect;
    Polynomial P;

    do
    {
        if (Incorrect.size() < 10)
        {
            for (int i = 0; i < Incorrect.size(); i++)
            {
                RndInterval I = Incorrect.at(i);
                X.push_back(I);
            }
        }
        else
        {
            for (int i = 0; i < 10; i++)
            {
                RndInterval I = Incorrect.at(floor(i * Incorrect.size() / 11));
                X.push_back(I);
            }
        }
        printf("Sample size: %ld\n", X.size());
        printf("Generating RndIntervals...\n");
        vector<RndInterval> L = CalcRndIntervals(X);

        printf("Generating RedIntervals...\n");
        vector<RndInterval> L2 = CalcRedIntervals(L);

        printf("Generating Polynomial...\n");
        P = GeneratePolynomial(L2);

        Incorrect = Verify(L2, P);
        if (Incorrect.size() > 0)
        {
            printf("FAILED IN SAMPLE: %ld\n", Incorrect.size());
            return 0;
        }
        Incorrect = Verify(Test, P);
    } while (Incorrect.size() > 0);

    for (int i = 0; i < P.termsize; i++)
    {
        printf("%5.60f\n", P.coefficients.at(i));
    }

    FILE *fptr = fopen("dump/poly.txt", "w");
    for (int i = 0; i < P.termsize; i++)
    {
        fprintf(fptr, "%5.60f\n", P.coefficients.at(i));
    }
    fclose(fptr);

    printf("Finished!\n");
}