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

vector<RndInterval> GenerateFloatSample(float start, float end, size_t cap, size_t step)
{
    FILE *fptr;
    fptr = fopen("dump/GenFloatSample.txt", "w");

    vector<RndInterval> X;
    size_t n = 0;

    for (float f = start; f < end;
         f = nextafterf(f, INFINITY))
    {
        if (n > cap)
        {
            break;
        }

        n++;

        RndInterval I;
        I.x_orig = f;
        I.x_rr = RangeReduction(f);
        X.push_back(I);

        fprintf(fptr, "%4.30f %4.30f\n", f, I.x_rr);

        for (int i = 0; i < step; i++)
        {
            f = nextafterf(f, INFINITY);
        }
    }
    printf("sample size = %ld\n", n);
    return X;
}

vector<RndInterval> CalcRndIntervals(vector<RndInterval> X)
{
    FILE *fptr;
    fptr = fopen("dump/CalcRndIntervals.txt", "w");

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

        mpfr_set_flt(x_mpfr, X.at(i).x_orig, MPFR_RNDN);
        mpfr_log2(y_mpfr, x_mpfr, MPFR_RNDN);

        y = mpfr_get_flt(y_mpfr, MPFR_RNDN);

        float lfloat = nextafterf(y, -INFINITY);
        float ufloat = nextafterf(y, +INFINITY);

        l = ((double)y + (double)lfloat) / 2.0;
        u = ((double)y + (double)ufloat) / 2.0;

        while ((float)l != (float)y)
        {
            l = nextafter(l, u);
        }
        l = nextafter(l, u);

        assert((float)l == (float)y);

        while ((float)u != (float)y)
        {
            u = nextafter(u, l);
        }
        u = nextafter(u, l);

        assert((float)u == (float)y);

        if (l > u)
        {
            printf("x = %4.15f y = %4.15f l = %4.15f u = %4.15f\n", X.at(i).x_orig, y, l, u);
            double temp = l;
            l = u;
            u = temp;
        }
        assert(l < u);

        // fprintf(fptr, "%.24e %.24e %.24e %.24e \n", x, y, l, h);
        fprintf(fptr, "%4.15f %4.15f %4.40f %4.40f \n", (double)X.at(i).x_orig, (double)y, l, u);
        I.y = y;
        I.l = l;
        I.u = u;

        L.push_back(I);
    }
    return L;
}

vector<RndInterval> CalcRedIntervals(vector<RndInterval> L)
{
    FILE *fptr;
    fptr = fopen("dump/CalcRedIntervals.txt", "w");
    float yp;
    double lp, up;

    mpfr_t yp_mpfr, xrr_mpfr;
    mpfr_inits2(200, yp_mpfr, xrr_mpfr, NULL);

    vector<RndInterval> L2;

    for (size_t i = 0; i < L.size(); i++)
    {
        RndInterval I;
        I.x_orig = L.at(i).x_orig;
        I.x_rr = L.at(i).x_rr;

        I.y = L.at(i).y;
        I.l = L.at(i).l;
        I.u = L.at(i).u;

        lp = InverseOutputCompensation(L.at(i).x_orig, L.at(i).l);
        up = InverseOutputCompensation(L.at(i).x_orig, L.at(i).u);

        // printf("lp = %4.15f up = %4.15f\n", lp, up);
        // printf("lp = %4.15f up = %4.15f\n", lp, up);
        // printf("yp = %4.15f\n", yp);
        while ((float)OutputCompensation(L.at(i).x_orig, lp) != L.at(i).y)

        // while ((float)OutputCompensation(L.at(i).x_orig, lp) != L.at(i).y)
        {
            for (int j = 0; j < 3; j++)
            {
                lp = nextafter(lp, up);
            }
            // lp = nextafter(lp, up);
            // lp += 0.0000001;
        }
        // while (OutputCompensation(L.at(i).x_orig, up) > L.at(i).u || OutputCompensation(L.at(i).x_orig, up) < L.at(i).l)
        while ((float)OutputCompensation(L.at(i).x_orig, up) != L.at(i).y)
        {
            for (int j = 0; j < 3; j++)
            {
                up = nextafter(up, lp);
            }
            // up = nextafter(up, lp);
            // up -= 0.0000001;
        }
        // printf("lp = %4.15f up = %4.15f\n", lp, up);
        // printf("lp = %4.15f up = %4.15f\n", lp, up);
        // if (lp >= up)
        // {
        //     printf("x = %4.15f lp = %4.15f up = %4.15f\n", L.at(i).x_orig, lp, up);
        //     double temp = lp;
        //     lp = up;
        //     up = temp;
        // }
        assert(lp < up);

        fprintf(fptr, "%4.15f %4.40f %4.40f %4.40f \n", (double)L.at(i).x_orig, L.at(i).x_rr, lp, up);
        I.lp = lp;
        I.up = up;

        L2.push_back(I);
    }
    return L2;
}

Polynomial GeneratePolynomial(vector<RndInterval> L2)
{

    // int termsize = 3;
    // for (int termsize = 29; termsize < 30; termsize++)
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

        /* then constraints one by one */
        for (int i = 0; i < L2.size(); i++)
        // for (int i = 1; i < 50; i++)
        {
            DSVectorRational row1(termsize);
            Rational acc(1.0);
            Rational xval(L2.at(i).x_rr);
            row1.add(0, 1.0);

            for (int j = 1; j < termsize; j++)
            {
                acc = acc * xval;
                row1.add(j, acc);
            }
            double lbnd = (L2.at(i).lp);
            double ubnd = (L2.at(i).up);
            // print x
            // mpfr_out_str(stdout, 10, 0, L.at(i).x, MPFR_RNDN);
            mysoplex.addRowRational(LPRowRational(lbnd, row1, ubnd));
        }

        mysoplex.writeFileRational("dump/f32_dump_rational.lp", NULL, NULL, NULL);
        mysoplex.writeFileReal("dump/f32_dump_real.lp", NULL, NULL, NULL);

        SPxSolver::Status stat;
        // cout << "Solving... " << termsize << "\n";

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
        // std::cout << "Status: " << stat << std::endl;
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

vector<RndInterval> VerifyConsistant(vector<RndInterval> L2, Polynomial P)
{
    size_t n = 0;
    size_t correct = 0;

    for (size_t i = 0; i < L2.size(); i++)
    {
        n++;
        float f = L2.at(i).x_orig;

        // double x = (double)f;
        double x_rr = RangeReduction(f);
        double eval = EvaulutePoly(P, x_rr);

        float y = OutputCompensation(f, eval);

        mpfr_t y_mpfr;
        mpfr_init2(y_mpfr, 200);
        mpfr_set_flt(y_mpfr, f, MPFR_RNDN);
        mpfr_log2(y_mpfr, y_mpfr, MPFR_RNDN);

        float oracle = (float)mpfr_get_d(y_mpfr, MPFR_RNDN);

        if (y == oracle)
        {
            correct++;
        }
        else
        {
            // L2.at(i).l
            // print x y oracle xrr eval y l u lp up
            // printf("failed: x = %4.30f y = %4.30f oracle = %4.30f\n", f, y, oracle);
            // printf("l = %4.30f u = %4.30f lp  = %4.30f up = %4.30f\n", L2.at(i).l, L2.at(i).u, L2.at(i).lp, L2.at(i).up);
            // printf("x_rr = %4.30f eval = %4.30f\n\n", x_rr, eval);
            // print("")
        }
    }

    printf("In sample: correct = %ld incorrect = %ld\n", correct, n - correct);
    return L2;
}
vector<RndInterval> VerifyContinuous(vector<RndInterval> L2, Polynomial P, float start, float end, size_t cap)
{
    size_t n = 0;
    size_t correct = 0;

    for (float f = start; f < end;
         f = nextafterf(f, INFINITY))
    {
        if (cap > 0 & n > cap)
        {
            break;
        }

        n++;

        double x = (double)f;
        double x_rr = RangeReduction(f);
        double eval = EvaulutePoly(P, x_rr);

        float y = OutputCompensation(f, eval);

        mpfr_t y_mpfr;
        mpfr_init2(y_mpfr, 200);
        mpfr_set_d(y_mpfr, x, MPFR_RNDN);
        mpfr_log2(y_mpfr, y_mpfr, MPFR_RNDN);

        float oracle = mpfr_get_flt(y_mpfr, MPFR_RNDN);

        if (y == oracle)
        {
            correct++;
        }
        else
        {
            // printf("continuous Failed: x = %4.30f y = %4.30f oracle = %4.30f x_rr = %4.30f eval = %4.30f\n", x, y, oracle, x_rr, eval);
            RndInterval I;
            I.x_orig = f;
            I.x_rr = x_rr;

            L2.push_back(I);
            return L2;
        }
    }
    printf("Full range: correct = %ld incorrect = %ld\n", correct, n - correct);
    return L2;
}

int main(int argc, char *argv[])
{
    // 0 : dont use files at all
    // 1 : send results to files
    // 2 : read results from files if available

    int save_to_file = 0;
    // check if first command line arg is -y
    if (argc == 2)
    {
        if (strcmp(argv[1], "1") == 0)
        {
            save_to_file = 1;
        }
    }

    printf("Generating FloatSample...\n");
    vector<RndInterval> X = GenerateFloatSample(SAMPLE_START, SAMPLE_END, 1000000, 1234);
    size_t last = X.size();

    Polynomial P;
    for (;;)
    {

        printf("Generating RndIntervals...\n");
        vector<RndInterval> L = CalcRndIntervals(X);

        printf("Generating RedIntervals...\n");
        vector<RndInterval> L2 = CalcRedIntervals(L);

        printf("Generating Polynomial...\n");
        P = GeneratePolynomial(L2);

        FILE *fptr;
        fptr = fopen("dump/poly.txt", "w");
        for (int i = 0; i < P.termsize; i++)
        {
            fprintf(fptr, "%5.60f\n", P.coefficients.at(i));
            // printf("%5.60f\n", P.coefficients.at(i));
        }

        X = VerifyConsistant(L2, P);
        X = VerifyContinuous(X, P, SAMPLE_START, SAMPLE_END, -1);
        if (X.size() == last)
        {
            break;
        }
        printf("X size = %ld\n", X.size());
        last = X.size();
    }
    X = VerifyContinuous(X, P, SAMPLE_START, SAMPLE_END, -1);

    printf("Finished!\n");
}