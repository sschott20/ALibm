#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <inttypes.h>
#include <math.h>
int bit_return(int a, int loc)
{
    int buf = a & 1 << loc;

    if (buf == 0)
        return 0;
    else
        return 1;
}

int print_float(float a)
{

    // 11000010111011010100000000000000
    //  1 sign bit | 8 exponent bit | 23 fraction bits
    int *b;
    b = (int *)&a;

    int i;
    for (i = 31; i >= 0; i--)
    {
        if (i == 30 || i == 22)
        {
            printf(" ");
        }
        printf("%d", bit_return(*b, i));
    }
    printf(" %.150f", a);
    printf("\n");
    return 0;
}

int main(void)
{

    FILE *fptr;
    fptr = fopen("floats.txt", "w");
    float current_float = 0;

    mpfr_t float_300;
    mpfr_init2(float_300, 300);

    for (int i = 0; i < 1'704; i++)
    // for (int i = 0; i < 947'912'704; i++)
    {

        mpfr_set_flt(float_300, current_float, MPFR_RNDN);
        mpfr_log2(float_300, float_300, MPFR_RNDN);

        int *b;
        b = (int *)&current_float;
        for (int i = 31; i >= 0; i--)
        {
            if (i == 30 || i == 22)
            {
                fprintf(fptr, " ");
            }
            fprintf(fptr, "%d", bit_return(*b, i));
        }
        fprintf(fptr, ", ");

        fprintf(fptr, "%.24e, ", current_float);

        float y = mpfr_get_flt(float_300, MPFR_RNDN);
        float l, h;
        l = nextafterf(y, 0.0f);
        h = nextafterf(y, 1.0f);

        fprintf(fptr, "%.24e ", y);
        fprintf(fptr, "%.24e ", l);
        fprintf(fptr, "%.24e ", h);
        fprintf(fptr, "\n");

        current_float = nextafterf(current_float, 1.0f);
    }
    mpfr_clear(float_300);
}
