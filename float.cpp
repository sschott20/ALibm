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
    float x;
    mpfr_t lower;
    mpfr_t upper;
};

struct RedInterval
{
    mpfr_t x;
    mpfr_t lower;
    mpfr_t upper;
    float x_orig;
    mpfr_t l_orig;
    mpfr_t u_orig;
};

struct Polynomial
{
    int termsize;
    vector<double> coefficients;
};

union float_bits {
    float f;
    uint32_t i;
};


vector<RndInterval> GenerateFloatSample()
{
    vector<RndInterval> X;
    size_t n = 0;

    for (float floatValue = 1; floatValue <= 1.1;
         floatValue = nextafterf(floatValue, INFINITY))
    {
        n++;
        RndInterval I;
        I.x = floatValue;
        X.push_back(I);
        for (int i = 0; i < 100; i++)
        {
            floatValue = nextafterf(floatValue, INFINITY);
        }

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

    mpfr_inits2(200, y_mpfr, x_mpfr, l_mpfr, h_mpfr, middle, NULL);

    vector<RndInterval> L;

    for (size_t i = 0; i < X.size(); i++)
    {
        x = X.at(i).x;

        RndInterval I;

        mpfr_init2(I.lower, 200);
        mpfr_init2(I.upper, 200);

        I.x = x;

        mpfr_set_flt(x_mpfr, x, MPFR_RNDN);
        mpfr_log2(y_mpfr, x_mpfr, MPFR_RNDN);

        y = mpfr_get_flt(y_mpfr, MPFR_RNDN);
        h = nextafterf(y, +INFINITY);
        l = nextafterf(y, -INFINITY);

        mpfr_set_flt(middle, y, MPFR_RNDN);
        mpfr_set_flt(l_mpfr, l, MPFR_RNDN);
        mpfr_set_flt(h_mpfr, h, MPFR_RNDN);

        mpfr_add(l_mpfr, l_mpfr, middle, MPFR_RNDN);
        mpfr_div_ui(l_mpfr, l_mpfr, 2, MPFR_RNDN);

        // should be while loop
        if (mpfr_get_flt(l_mpfr, MPFR_RNDN) != y)
        {
            mpfr_nextabove(l_mpfr);
        }
        mpfr_set(I.lower, l_mpfr, MPFR_RNDN);
        assert(mpfr_get_flt(I.lower, MPFR_RNDN) == y);

        mpfr_add(h_mpfr, h_mpfr, middle, MPFR_RNDN);
        mpfr_div_ui(h_mpfr, h_mpfr, 2, MPFR_RNDN);

        // should be while loop
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

    for (size_t i = 1; i < L.size(); i++)
    {
        // fprintf(fptr, "%4.30Lf %4.30Lf\n", mpfr_get_ld(L.at(i).upper, MPFR_RNDN), mpfr_get_ld(L.at(i).lower, MPFR_RNDN));
        // mpfr_out_str(fptr, 10, 0, L.at(i).lower, MPFR_RNDD);
        mpfr_t x_rr, int_part, exp_mpfr, x2_mpfr;
        mpfr_t alpha, beta;
        mpfr_inits2(200, x_rr, int_part, exp_mpfr, x2_mpfr, NULL);
        mpfr_inits2(200, alpha, beta, NULL);

        mpfr_log2(exp_mpfr, x2_mpfr, MPFR_RNDN);
        mpfr_floor(exp_mpfr, exp_mpfr);
        long exp = mpfr_get_si(exp_mpfr, MPFR_RNDZ);
        exp = exp + 2;
        // fprintf(fptr, "%ld\n", exp);
        mpfr_set_flt(exp_mpfr, 2, MPFR_RNDN);
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
        mpfr_inits2(200, I.lower, I.upper, I.x, I.l_orig, I.u_orig, NULL);

        mpfr_set(I.x, x_rr, MPFR_RNDN);
        mpfr_set(I.lower, alpha, MPFR_RNDN);
        mpfr_set(I.upper, beta, MPFR_RNDN);

        I.x_orig = L.at(i).x;
        mpfr_set(I.l_orig, L.at(i).lower, MPFR_RNDN);
        mpfr_set(I.u_orig, L.at(i).upper, MPFR_RNDN);

        fprintf(fptr, "%4.40f \n", L.at(i).x);
        mpfr_fprintf(fptr, "%3.48Rf %3.48Rf %3.48Rf \n", I.x, I.lower, I.upper);

        fprintf(fptr, "\n");
        L2.push_back(I);
    }
    return L2;
}

Polynomial GeneratePolynomial(vector<RedInterval> L)
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
            Rational xval(mpfr_get_d1(L.at(i).x));
            row1.add(0, 1.0);

            for (int j = 1; j < termsize; j++)
            {
                acc = acc * xval;
                row1.add(j, acc);
            }
            double lbnd = mpfr_get_d1(L.at(i).lower);
            double ubnd = mpfr_get_d1(L.at(i).upper);
            // print x
            // mpfr_out_str(stdout, 10, 0, L.at(i).x, MPFR_RNDN);
            mysoplex.addRowRational(LPRowRational(lbnd, row1, ubnd));
        }

        mysoplex.writeFileRational("dump_rational.lp", NULL, NULL, NULL);
        mysoplex.writeFileReal("dump_real.lp", NULL, NULL, NULL);

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
double EvaulutePoly(Polynomial P, float x)
{
    double xval = (double)x;
    double acc = 0.0;
    double power = 1.0;
    for (int i = 0; i < P.termsize; i++)
    {
        // printf("%ld", P.coefficients.at(i));
        acc = acc + P.coefficients.at(i) * power;
        power = power * x;
    }
    return acc;
}

void Verify(vector<RedInterval> L2, Polynomial P)
{
    vector<RndInterval> X;
    size_t n = 0;
    size_t correct = 0;
    for (float floatValue = 1; floatValue <= 1.1; floatValue = nextafterf(floatValue, INFINITY))
    {
        n++;
        float range_reduced = floatValue;
        float eval = (float)EvaulutePoly(P, floatValue);
        mpfr_t y;
        mpfr_inits2(200, y, NULL);
        mpfr_set_flt(y, floatValue, MPFR_RNDN);
        mpfr_log2(y, y, MPFR_RNDN);

        float oracle = mpfr_get_flt(y, MPFR_RNDN);
        printf("Float: %6.60f \n ", eval - floatValue);
        if (eval != oracle)
        {

            // printf("Float: %6.60f Eval: %6.60f Oracle: %6.50f\n ", floatValue, eval, oracle);
        }
        else
        {
            correct += 1;
        }

        if (n > 1000)
        {
            break;
        }

        // cout << "Float: " << floatValue << " Eval: " << eval << "\n";
    }
    printf("Correct: %ld Total: %ld\n", correct, n);
    return;
}
int main()
{

    printf("Generating FloatSample...\n");
    vector<RndInterval> X = GenerateFloatSample();

    printf("Generating RndIntervals...\n");
    vector<RndInterval> L = CalcRndIntervals(X);

    printf("Generating RedIntervals...\n");
    vector<RedInterval> L2 = CalcRedIntervals(L);

    printf("Generating Polynomial...\n");

    Polynomial P = GeneratePolynomial(L2);
    for (int i = 0; i < P.termsize; i++)
    {
        printf("%5.60f\n", P.coefficients.at(i));
    }

    // Verify(L2, P);

    printf("Finished!\n");
}