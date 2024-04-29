#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;

struct Polynomial
{
    int termsize;
    vector<double> coefficients;
};
typedef union
{
    float f;
    unsigned x;
} floatX;

typedef union
{
    double d;
    unsigned long x;
} doubleX;

unsigned long flog2_coef[5][7] = {

    {
        13840088983442381949U,
        4625347124714385442U,
        13854123740975661436U,
        4634414482658525937U,
        13858134739117318847U,
        4631380976507741609U,
        13846118394287531079U,

    },
    {
        13839950078055808653U,
        4624734434826564558U,
        13853250229213416370U,
        4632490643537813669U,
        13855758816733519427U,
        4628755719687040412U,
        13842904545747925551U,

    },
    {
        13839813994467456793U,
        4624047845104605039U,
        13851984055115512867U,
        4630870402802447671U,
        13853712686755463775U,
        4625998390763997241U,
        13839806877980237190U,

    },
    {
        13839685728680924647U,
        4623451467538231135U,
        13850828237194977420U,
        4629650078793719008U,
        13851571232100972718U,
        4623345689582687453U,
        13836476175883689232U,

    },
    {
        13839589662576611435U,
        4623033249594088492U,
        13850069656578378951U,
        4628182459009314188U,
        13849974111505167414U,
        4621491791609445057U,
        13834307695743617449U,

    }};
Polynomial readPolynomialFromFile(const std::string &filename);
double EvaulutePoly(int index, double xval);
Polynomial InitPoly(int index);