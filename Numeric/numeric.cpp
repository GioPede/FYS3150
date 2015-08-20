//various libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[])
{
    FILE * data;

    data = fopen("data.txt", "w+");

    int N = atoi(argv[1]);
    double dsum_up, dsum_down, M;
    float fsum_up, fsum_down;
    for(int j = 0; j < N; j++)
    {
        M = pow(10,j);
        dsum_up = 0;
        dsum_down = 0;
        fsum_up = 0;
        fsum_down = 0;

        for (double i = 1; i < M; i++)
        {
            dsum_up += 1.0 / i;
            fsum_up += 1.0 / i;
            dsum_down += 1.0 / (M-i);
            fsum_down += 1.0 / (M-i);
        }
        fprintf(data, "%g,  %.30g, %.30g\n", M, (fsum_down-fsum_up > 0 ? fsum_down-fsum_up : fsum_up - fsum_down),
                                                (dsum_down-dsum_up > 0 ? dsum_down-dsum_up : dsum_up - dsum_down));
        printf("%d\n",j);
        /*
        printf("Float Up: %.20g\n", fsum_up);
        printf("Float Down: %.20g\n", fsum_down);
        printf("Difference: %.20g\n\n", fsum_up-fsum_down);

        printf("Double Up: %.20g\n", dsum_up);
        printf("Double Down: %.20g\n", dsum_down);
        printf("Difference: %.20g\n\n", dsum_up-dsum_down);
        */
    }
    fclose(data);
}
