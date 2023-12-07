#include <stdio.h>
#include <iostream>
#include <inttypes.h>
#include <math.h>
#include <random>
#include <cstdlib>
#include <vector>

using namespace std;

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


int main()
{
    vector<RndInterval> X;
    for (int i = 0; i < 123; i++)
    {
        RndInterval x1 = {x_orig : i, x_rr : 2};
        X.push_back(x1);
    }

    write_rnd_interval_to_file(X, "iotest.bin");
    vector<RndInterval> X2 = read_rnd_interval_from_file("iotest.bin");

    for (int i = 0; i < X2.size(); i++)
    {
        printf("All fields of X: %f %f %f %f %f %f %f %f\n", X.at(i).x_orig, X.at(i).y, X.at(i).l, X.at(i).u, X.at(i).x_rr, X.at(i).yp, X.at(i).lp, X.at(i).up);
    }
    return 0;
}