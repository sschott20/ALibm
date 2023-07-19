#include <stdio.h>

#include <gmp.h>
#include <mpfr.h>
#include <math.h>
int main(void)
{
    float i = .0000001123f;
    float *g = &i;
    mpfr_t t;

    mpfr_init2(t, 200);
    mpfr_set_d(t, *g, MPFR_RNDD);
    mpfr_log2(t, t, MPFR_RNDD);
    float result = mpfr_get_flt(t, MPFR_RNDN);
    printf("%f\n", result);

    mpfr_clear(t);
    mpfr_free_cache();

    // float a = 0.0f;
    // printf("%e\n", a);
    // float b = nextafterf(a, 1.0f);
    // printf("%e\n", b);
    return 0;
}