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
// rr to [0.5, 1)
double RangeReduction(half x)
{
    // reduces x to [0.5, 1)
    int exp;
    double sig = half_float::frexp(x, &exp);
    // exp -= 1;
    // sig *= 2;
    // print_binary_half(x);
    // std::cout << x << " makes " << sig << " * 2^" << exp - 1 << '\n';

    // return sig * 2;
    // printf("x = %4.15f sig = %4.15f exp = %d\n", (double)x, sig, exp);

    return sig;
}

double OutputCompensation(half x, double yp)
{
    int exp;
    double sig = half_float::frexp(x, &exp);
    // sig *= 2;
    // exp -= 1;
    // return yp + exp - 1;
    // printf("x = %4.15f yp = %4.15f, exp = %d\n", (double)x, yp, exp);

    return yp + exp;
}

double InverseOutputCompensation(half x, double yp)
{
    int exp;
    double sig = half_float::frexp(x, &exp);
    // sig *= 2;
    // exp -= 1;
    // return yp - (exp + 1);
    // printf("x = %4.15f yp = %4.15f, exp = %d ret = %4.15f\n", (double)x, yp, exp, yp - exp + 1);
    return yp - exp;
}

half EvaluateFunction(mpfr_t y, double x)
{
    mpfr_set_d(y, x, MPFR_RNDN);
    mpfr_log2(y, y, MPFR_RNDN);
    half h = (half)mpfr_get_d(y, MPFR_RNDN);
    return h;
}

vector<RndInterval> GenerateFloatSample()
{
    vector<RndInterval> X;
    size_t n = 0;
    half h;
    half zero(0.0);
    half inf((0x7BFF));
    for (h = half_float::nextafter(zero, inf); h < 1; h = half_float::nextafter(h, inf))
    {
        n++;
        RndInterval I;
        I.x_orig = h;
        I.x_rr = RangeReduction(h);
        X.push_back(I);
        for (int i = 0; i < 100; i++)
        {
            h = half_float::nextafter(h, inf);
        }
    }
    printf("sample size = %ld\n", n);
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

        mpfr_set_flt(x_mpfr, X.at(i).x_orig, MPFR_RNDN);
        mpfr_log2(y_mpfr, x_mpfr, MPFR_RNDN);

        y = (half)mpfr_get_d(y_mpfr, MPFR_RNDN);

        half inf(half_float::half(0x7BFF));

        half lhalf(half_float::nextafter((half)y, -inf));
        half uhalf(half_float::nextafter((half)y, inf));

        l = (y + (double)lhalf) / 2;
        u = (y + (double)uhalf) / 2;

        while ((half)l != (half)y)
        {
            // l += 0.00000001;
            // l = half_float::nextafter(l, inf);
            // for (int i = 0; i < 10; i++)
            // {
            l = nextafterf(l, INFINITY);
            // }
        }
        assert((half)l == (half)y);

        while ((half)u != (half)y)
        {
            // u -= 0.00000001;
            // u = half_float::nextafter(u, -inf);
            // for (int i = 0; i < 10; i++)
            // {
            u = nextafterf(u, -INFINITY);
            // }
        }
        assert((half)u == (half)y);
        // if (l > u){
        //     // print l u y x
        //     printf("l = %4.15f u = %4.15f y = %4.15f x = %4.15f\n", l, u, (double) y, (double)X.at(i).x_orig);
        // }
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

    // mpfr_t y_mpfr, x_mpfr, l_mpfr, h_mpfr, middle;
    // mpfr_inits2(200, y_mpfr, x_mpfr, l_mpfr, h_mpfr, middle, NULL);
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
            // lp += 0.00000001;
            lp = nextafterf(lp, INFINITY);
        }

        while (OutputCompensation(X.at(i).x_orig, up) > X.at(i).u)
        {
            // up -= 0.00000001;
            up = nextafterf(up, -INFINITY);
        }

        assert(lp <= up);
        fprintf(fptr, "%4.15f %4.40f %4.40f %4.40f \n", (double)X.at(i).x_orig, X.at(i).x_rr, lp, up);
        I.lp = lp;
        I.up = up;

        L.push_back(I);
        // printf("\n");
    }
    return L;
}

Polynomial GeneratePolynomial(vector<RndInterval> L)
{

    // int termsize = 3;
    // for (int termsize = 29; termsize < 30; termsize++)
    for (int termsize = 7; termsize < 30; termsize++)
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

        /* then constraints one by one */
        for (int i = 0; i < L.size(); i++)
        // for (int i = 1; i < 50; i++)
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
            // print x
            // mpfr_out_str(stdout, 10, 0, L.at(i).x, MPFR_RNDN);
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
        // output stat
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

void Verify(vector<RndInterval> L2, Polynomial P)
{
    size_t n = 0;
    size_t correct = 0;
    half inf(half_float::half(0x7BFF));

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

        // mpfr_set_d(y_mpfr, x, MPFR_RNDN);
        // mpfr_log2(y_mpfr, y_mpfr, MPFR_RNDN);

        // half oracle = (half)mpfr_get_d(y_mpfr, MPFR_RNDN);
        half oracle = EvaluateFunction(y_mpfr, x);

        if (y != oracle)
        {
            printf("Failed in sample! %4.10f %4.10f, %4.10f\n", (double)h, (double)oracle, (double)y);
        }
        else
        {
            // printf("Correct!\n");
            correct += 1;
        }

        // if (n > 10000)
        // {
        //     break;
        // }

        // cout << "Float: " << floatValue << " Eval: " << eval << "\n";
    }

    printf("In sample, Correct: %ld Incorrect: %ld\n", correct, n - correct);
    return;
}

vector<RndInterval> VerifyAll(vector<RndInterval> L2, Polynomial P)
{
    size_t n = 0;
    size_t correct = 0;
    half inf(half_float::half(0x7BFF));
    // half h(half_float::half(0xFC00));
    half h(half_float::half(0x0000));

    vector<RndInterval> Incorrect;
    for (;;)
    {
        h = half_float::nextafter(h, inf);
        if (h == inf)
        {
            break;
        }
        double x = (double)h;
        double range_reduced = RangeReduction(h);
        double eval = EvaulutePoly(P, range_reduced);
        half y(OutputCompensation(h, eval));

        mpfr_t y_mpfr;
        mpfr_inits2(200, y_mpfr, NULL);
        // mpfr_set_d(y_mpfr, x, MPFR_RNDN);
        // mpfr_log2(y_mpfr, y_mpfr, MPFR_RNDN);

        // half oracle = (half)mpfr_get_d(y_mpfr, MPFR_RNDN);
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
        n++;
    }

    printf("Out of sample, Correct: %ld Incorrect: %ld\n", correct, n - correct);
    return Incorrect;
}

int main()
{

    // half x(100.5);
    // double xrr = RangeReduction(x);
    // printf("x = %4.15f, xrr = %4.15f\n", (double)x, xrr);

    printf("Generating FloatSample...\n");
    vector<RndInterval> X = GenerateFloatSample();

    vector<RndInterval> Incorrect;
    Polynomial P;

    do
    {
        // new sample adds 10 more values to X from incorrect
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

        Verify(L2, P);
        Incorrect = VerifyAll(L2, P);
    } while (Incorrect.size() > 0);

    for (int i = 0; i < P.termsize; i++)
    {
        printf("%5.60f\n", P.coefficients.at(i));
    }

    printf("Finished!\n");
}