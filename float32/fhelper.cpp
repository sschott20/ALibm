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
#include "fhelper.hpp"

using namespace std;
using namespace soplex;

vector<RndInterval> GenerateFloatSample(int sample_size, float min, float max)
{
    vector<RndInterval> X;
    size_t n = 0;
    float h;

    if (sample_size == -1)
    {
        for (h = min; h < max; h = nextafterf(h, max))
        {
            if (h == INFINITY || h == -INFINITY)
            {
                break;
            }
            if (ComputeSpecialCase(h) == -1)
            {
                continue;
            }
            n++;
            RndInterval I;
            I.x_orig = h;
            I.x_rr = RangeReduction(h);
            X.push_back(I);
            // printf("Sample size: %ld\n", n);

            for (int i = 0; i < 6; i++)
            {
                h = nextafterf(h, max);
            }
        }
    }
    else
    {
        // long skip = (1 << 16) / sample_size;
        unsigned long long skip = (1 << 31) / sample_size;

        for (h = min; h < max; h = nextafterf(h, max))
        {
            if (h == INFINITY || h == -INFINITY)
            {
                break;
            }
            if (ComputeSpecialCase(h) == -1)
            {
                continue;
            }

            n++;
            RndInterval I;
            // float rand_h = h + (float)rand() / (float)(RAND_MAX / (skip));
            I.x_orig = h;
            I.x_rr = RangeReduction(h);
            X.push_back(I);

            // h += skip;
            for (int i = 0; i < skip; i++)
            {
                h = nextafterf(h, max);
            }
        }
    }
    printf("Sample size: %ld\n", X.size());
    return X;
}

vector<RndInterval> CalcRndIntervals(vector<RndInterval> X)
{
    FILE *fptr;
    fptr = fopen("../dump/FRndInterval.txt", "w");

    float y;
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

        float lhalf(nextafterf((float)y, -INFINITY));
        float uhalf(nextafterf((float)y, INFINITY));

        l = ((double)y + (double)lhalf) / 2;
        u = ((double)y + (double)uhalf) / 2;

        // printf("x: %4.40f y: %4.40f l: %4.40f u: %4.40f\n", (double)X.at(i).x_orig, (double)y, l, u);
        while ((float)l != (float)y)
        {
            l = nextafter(l, INFINITY);
        }
        // l = nextafter(l, INFINITY);

        assert((float)l == y);

        while ((float)u != y)
        {
            u = nextafter(u, -INFINITY);
        }
        // u = nextafter(u, -INFINITY);

        assert((float)u == (float)y);
        assert(l < u);

        // fprintf(fptr, "%4.40f %4.40f %4.40f %4.40f \n", (double)X.at(i).x_orig, (double)y, l, u);
        I.y = y;
        I.l = l;
        I.u = u;

        L.push_back(I);
    }
    std::sort(L.begin(), L.end());
    L.erase(std::unique(L.begin(), L.end()), L.end());
    return L;
}

vector<RndInterval> CalcRedIntervals(vector<RndInterval> X)
{
    FILE *fptr;
    fptr = fopen("../dump/FRedInterval.txt", "w");
    float yp;
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

        // lp = InverseOutputCompensation(X.at(i).x_orig, X.at(i).l);
        // up = InverseOutputCompensation(X.at(i).x_orig, X.at(i).u);
        // while (OutputCompensation(X.at(i).x_orig, lp) < X.at(i).l)
        // {
        //     lp = nextafter(lp, INFINITY);
        // }

        // while (OutputCompensation(X.at(i).x_orig, up) > X.at(i).u)
        // {
        //     up = nextafter(up, -INFINITY);
        // }

        GuessInitial(X.at(i).x_orig, lp, up);
        unsigned long long step = 0x8000000000000llu;
        double oc;
        // while (step > 0)
        while (step >= 0x0000000000800llu)

        {
            doubleX dx;
            dx.d = lp;
            if (dx.d >= 0)
            {
                dx.x -= step;
            }
            else
            {
                dx.x += step;
            }

            oc = OutputCompensation(I.x_orig, dx.d);
            if (oc >= I.l && oc <= I.u)
            {
                lp = dx.d;
            }
            else
            {
                step /= 2;
            }
        }
        step = 0x8000000000000llu;
        // while (step > 0)
        while (step >= 0x0000000000800llu)
        {
            doubleX dx;
            dx.d = up;
            if (dx.d >= 0)
            {
                dx.x += step;
            }
            else
            {
                dx.x -= step;
            }

            oc = OutputCompensation(I.x_orig, dx.d);
            if (oc > I.l && oc < I.u)
            {
                up = dx.d;
            }
            else
            {
                step /= 2;
            }
        }

        // printf("x: %4.40f xrr: %4.40f y: %4.40f l: %4.40f u: %4.40f lp: %4.40f up: %4.40f\n", I.x_orig, I.x_rr, I.y, I.l, I.u, lp, up);

        oc = OutputCompensation(I.x_orig, lp);
        // printf("oc: %4.40f\n", oc);
        assert((float)oc == I.y);
        oc = OutputCompensation(I.x_orig, up);
        // printf("oc: %4.40f\n", oc);
        assert((float)oc == I.y);

        assert(lp <= up);
        // fprintf(fptr, "%4.40f %4.40f %4.40f %4.40f %4.40f \n", (double)X.at(i).x_orig, X.at(i).x_rr, X.at(i).y, X.at(i).l, X.at(i).u);
        // fprintf(fptr, "%4.40f %4.40f %4.40f  \n", (double)X.at(i).yp, X.at(i).lp, X.at(i).up);
        I.lp = lp;
        I.up = up;

        L.push_back(I);
    }
    return L;
}

Polynomial GeneratePolynomial(vector<RndInterval> L)
{

    cout << "Solving... [1, 50] \n";
    for (int termsize = 1; termsize < 50; termsize++)
    {
        // printf("Termsize: %d\n", termsize);
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

        mysoplex.writeFileRational("../dump/FSoplexRational.lp", NULL, NULL, NULL);
        mysoplex.writeFileReal("../dump/FSoplexReal.lp", NULL, NULL, NULL);

        SPxSolver::Status stat;

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
                // P.coefficients.push_back(prim[i].doubleValue())
                // P.coefficients.push_back(mpq_get_d(*(prim[i].getMpqPtr())));
            }
            std::cout << "Status: " << stat << " with " << termsize << " terms" << std::endl;

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
        // std::cout << "Status: " << stat << std::endl;
    }

    Polynomial N;
    return N;
}

long double EvaulutePoly(Polynomial P, double xval)
{
    long double acc = 0.0;
    long double power = 1.0;

    for (int i = 0; i < P.termsize; i++)
    {

        acc += P.coefficients.at(i) * power;
        power *= xval;
    }

    return acc;
}

vector<RndInterval> Verify(vector<RndInterval> L2, Polynomial P, int debug = 0)
{
    size_t n = 0;
    size_t correct = 0;

    vector<RndInterval> Incorrect;
    for (size_t i = 0; i < L2.size(); i++)
    {
        n++;
        float h = L2.at(i).x_orig;

        double x = (double)h;
        double range_reduced = RangeReduction(h);
        long double eval = EvaulutePoly(P, range_reduced);
        double oc = OutputCompensation(h, eval);
        float y = (float)oc;

        mpfr_t y_mpfr;
        mpfr_inits2(200, y_mpfr, NULL);
        float oracle = EvaluateFunction(y_mpfr, x);
        if (y != oracle)
        {
            RndInterval I;
            I.x_orig = h;
            I.x_rr = range_reduced;
            Incorrect.push_back(I);
            if (debug)
            {
                printf("x: %4.40f y: %4.40f oracle: %4.40f\n", h, y, oracle);
                // printf("xrr: %4.40f eval: %4.40f\n", range_reduced, eval);
                // printf("lp: %4.40f up: %4.40f\n", L2.at(i).lp, L2.at(i).up);
            }
        }
        else
        {
            correct += 1;
        }
    }

    printf("In sample, Correct: %ld Incorrect: %ld\n", correct, n - correct);
    return Incorrect;
}

void print_poly(Polynomial P)
{
    FILE *fptr = fopen("../dump/FPoly.txt", "w");
    for (int i = 0; i < P.termsize; i++)
    {
        printf("%5.60f\n", P.coefficients.at(i));
        fprintf(fptr, "%5.60f\n", P.coefficients.at(i));
    }
    fclose(fptr);
}

void FullTest(Polynomial P, float min, float max)
{
    size_t n = 0;
    size_t correct = 0;

    vector<RndInterval> Incorrect;

    for (float h = min; h < max; h = nextafterf(h, max))
    {
        // printf("x: %4.70f\n", h);
        if (h == INFINITY || h == -INFINITY)
        {
            break;
        }
        if (ComputeSpecialCase(h) == -1)
        {
            continue;
        }
        n++;

        double x = (double)h;
        double range_reduced = RangeReduction(h);
        long double eval = EvaulutePoly(P, range_reduced);
        double oc = OutputCompensation(h, eval);
        float y = (float)oc;

        mpfr_t y_mpfr;
        mpfr_inits2(200, y_mpfr, NULL);
        float oracle = EvaluateFunction(y_mpfr, x);
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
        for (int i = 0; i < 5; i++)
        {
            h = nextafterf(h, max);
        }
    }

    printf("Full Test Results, Correct: %ld Incorrect: %ld\n", correct, n - correct);
}
