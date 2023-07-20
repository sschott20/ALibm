#include <stdio.h>
#include <iostream>
#include <gmp.h>
#include <mpfr.h>
#include <inttypes.h>
#include <math.h>
#include <bits/stdc++.h>
#include <random>
#include <soplex.h>

using namespace std;
using namespace soplex;

struct RndInterval
{
    float x;
    mpfr_t lower;
    mpfr_t upper;
};

struct RedInterval
{
    mpfr_t x;
    mpfr_t lower;
    mpfr_t upper;
};

struct Polynomial
{
    vector<mpfr_t> coefficients;
};
vector<RndInterval> GenerateFloatSample()
{
    vector<RndInterval> X;
    size_t n = 0;

    for (float floatValue = 100; floatValue <= 100.5; floatValue = nextafterf(floatValue, INFINITY))
    {
        n++;
        // printf("floatValue = %f\n", floatValue);
        RndInterval I;
        I.x = floatValue;
        X.push_back(I);
    }
    printf("sample size = %ld\n", n);
    return X;
}
vector<RndInterval> CalcRndIntervals(vector<RndInterval> X)
{
    FILE *fptr;
    fptr = fopen("CalcRndIntervals.txt", "w");

    float x = 0;
    float y, l, h;

    mpfr_t y_mpfr, x_mpfr, l_mpfr, h_mpfr, middle;

    mpfr_inits2(500, y_mpfr, x_mpfr, l_mpfr, h_mpfr, middle, NULL);

    vector<RndInterval> L;

    for (size_t i = 0; i < X.size(); i++)
    {
        // printf("%d\n", i);

        x = X.at(i).x;

        RndInterval I;

        mpfr_init2(I.lower, 500);
        mpfr_init2(I.upper, 500);

        I.x = x;

        mpfr_set_flt(x_mpfr, x, MPFR_RNDN);
        mpfr_log10(y_mpfr, x_mpfr, MPFR_RNDN);

        y = mpfr_get_flt(y_mpfr, MPFR_RNDN);
        h = nextafterf(y, +INFINITY);
        l = nextafterf(y, -INFINITY);

        mpfr_set_flt(middle, y, MPFR_RNDN);
        mpfr_set_flt(l_mpfr, l, MPFR_RNDN);
        mpfr_set_flt(h_mpfr, h, MPFR_RNDN);

        mpfr_add(l_mpfr, l_mpfr, middle, MPFR_RNDN);
        mpfr_div_ui(l_mpfr, l_mpfr, 2, MPFR_RNDN);

        if (mpfr_get_flt(l_mpfr, MPFR_RNDN) != y)
        {
            mpfr_nextabove(l_mpfr);
        }
        mpfr_set(I.lower, l_mpfr, MPFR_RNDN);
        assert(mpfr_get_flt(I.lower, MPFR_RNDN) == y);

        mpfr_add(h_mpfr, h_mpfr, middle, MPFR_RNDN);
        mpfr_div_ui(h_mpfr, h_mpfr, 2, MPFR_RNDN);
        if (mpfr_get_flt(h_mpfr, MPFR_RNDN) != y)
        {
            mpfr_nextbelow(h_mpfr);
        }
        mpfr_set(I.upper, h_mpfr, MPFR_RNDN);
        assert(mpfr_get_flt(I.upper, MPFR_RNDN) == y);

        // fprintf(fptr, "%.24e %.24e %.24e %.24e \n", x, y, l, h);
        fprintf(fptr, "%4.30f %4.30f %4.30Lf %4.30Lf\n", x, y, mpfr_get_ld(I.upper, MPFR_RNDN), mpfr_get_ld(I.lower, MPFR_RNDN));
        L.push_back(I);
    }
    return L;
}

vector<RedInterval> CalcRedIntervals(vector<RndInterval> L)
{
    FILE *fptr;
    fptr = fopen("CalcRedIntervals.txt", "w");

    vector<RedInterval> L2;

    for (size_t i = 0; i < L.size(); i++)
    {
        // fprintf(fptr, "%4.30Lf %4.30Lf\n", mpfr_get_ld(L.at(i).upper, MPFR_RNDN), mpfr_get_ld(L.at(i).lower, MPFR_RNDN));
        // mpfr_out_str(fptr, 10, 0, L.at(i).lower, MPFR_RNDD);
        mpfr_t x_rr, int_part, exp_mpfr, x2_mpfr;
        mpfr_t alpha, beta;
        mpfr_inits2(500, x_rr, int_part, exp_mpfr, x2_mpfr, NULL);
        mpfr_inits2(500, alpha, beta, NULL);

        mpfr_log10(exp_mpfr, x2_mpfr, MPFR_RNDN);
        mpfr_floor(exp_mpfr, exp_mpfr);
        long exp = mpfr_get_si(exp_mpfr, MPFR_RNDZ);
        exp = exp + 3;
        // fprintf(fptr, "%ld\n", exp);
        mpfr_set_flt(exp_mpfr, 10, MPFR_RNDN);
        mpfr_pow_si(exp_mpfr, exp_mpfr, exp, MPFR_RNDN);

        mpfr_set_flt(x_rr, L.at(i).x, MPFR_RNDN);
        mpfr_div(x_rr, x_rr, exp_mpfr, MPFR_RNDN);

        mpfr_sub_si(alpha, L.at(i).lower, exp, MPFR_RNDN);
        mpfr_sub_si(beta, L.at(i).upper, exp, MPFR_RNDN);

        assert(mpfr_cmp(alpha, beta) < 0);

        while (true)
        {
            // check if alpha minus exp is above L.at(i).lower
            mpfr_set_ld(exp_mpfr, exp, MPFR_RNDN);
            mpfr_add(exp_mpfr, alpha, exp_mpfr, MPFR_RNDN);
            if (mpfr_cmp(exp_mpfr, L.at(i).lower) > 0)
            {
                mpfr_nextabove(alpha);
            }
            else
            {
                break;
            }
        }
        while (true)
        {
            // check if alpha minus exp is above L.at(i).lower
            mpfr_set_ld(exp_mpfr, exp, MPFR_RNDN);
            mpfr_add(exp_mpfr, beta, exp_mpfr, MPFR_RNDN);
            if (mpfr_cmp(exp_mpfr, L.at(i).upper) < 0)
            {
                mpfr_nextabove(beta);
            }
            else
            {
                break;
            }
        }

        RedInterval I;
        mpfr_inits2(500, I.lower, I.upper, I.x, NULL);

        mpfr_set(I.x, x_rr, MPFR_RNDN);
        mpfr_set(I.lower, alpha, MPFR_RNDN);
        mpfr_set(I.upper, beta, MPFR_RNDN);

        fprintf(fptr, "%4.40f \n", L.at(i).x);
        mpfr_fprintf(fptr, "%3.48Rf %3.48Rf %3.48Rf \n", I.x, I.lower, I.upper);

        fprintf(fptr, "\n");
        L2.push_back(I);
    }
    return L2;
}

Polynomial GeneratePolynomial(vector<RedInterval> L)
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

    /* we first add variables */
    DSVectorRational dummycol(0);
    for (int i = 0; i < 3; i++)
    {
        auto column = LPColRational(1.0, dummycol, infinity, -infinity);
        mysoplex.addColRational(column);
        // auto column = LPColReal(1.0, dummycol, infinity, -infinity);
        // mysoplex.addColReal(column);
    }

    /* then constraints one by one */
    // for (int i = 0; i < L.size(); i++)
    for (int i = 0; i < 10; i++)

    {
        DSVector row1(2);
        for (int i = 0; i < 3; i++)
        {
            // calcuate power of x using L.at(i).x

            row1.add(0, 1.0);
        }

        mysoplex.addRowRational(LPRowRational(L.at(i).lower, row1, L.at(i).upper));
        // mysoplex.addRowRational(LPRow(mpfr_get_d1(L.at(i).lower), row1, mpfr_get_d1(L.at(i).upper)));
    }

    // DSVector row1(2);
    // row1.add(0, 1.0);
    // row1.add(1, 5.0);
    // row1.add(2, 5.0);
    // mysoplex.addRowReal(LPRow(100.0, row1, infinity));

    /* write LP in .lp format */
    mysoplex.writeFileReal("dump.lp", NULL, NULL, NULL);

    /* solve LP */
    SPxSolver::Status stat;
    DVector prim(2);
    DVector dual(1);
    stat = mysoplex.optimize();

    /* get solution */
    if (stat == SPxSolver::OPTIMAL)
    {
        mysoplex.getPrimalReal(prim);
        mysoplex.getDualReal(dual);
        std::cout << "LP solved to optimality.\n";
        std::cout << "Objective value is " << mysoplex.objValueReal() << ".\n";
        std::cout << "Primal solution is [" << prim[0] << ", " << prim[1] << "].\n";
        std::cout << "Dual solution is [" << dual[0] << "].\n";
    }
    Polynomial P;

    return P;
}

int main()
{

    printf("Generating FloatSample...\n");
    vector<RndInterval> X = GenerateFloatSample();

    printf("Generating RndIntervals...\n");
    vector<RndInterval> L = CalcRndIntervals(X);

    printf("Generating RedIntervals...\n");
    vector<RedInterval> L2 = CalcRedIntervals(L);

    printf("Generating Polynumial...\n");
    Polynomial P = GeneratePolynomial(L2);
}