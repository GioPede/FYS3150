//various libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    double dsum_up = 0, dsum_down = 0;
    float fsum_up = 0, fsum_down = 0;

    double N = 1000000;

    for (int i = 0; i < N; i++)
    {
        dsum_up += 1.0 / i;
        fsum_up += 1.0 / i;
        dsum_down += 1.0 / (N-i);
        fsum_down += 1.0 / (N-i);
    }

    printf("Float Up: %g\n", fsum_up);
    printf("Double Up: %g\n", dsum_up);
    printf("Float Up: %g\n", fsum_down);
    printf("Double Up: %g\n", dsum_down);
}
