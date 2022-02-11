#include <stdio.h>

int main()
{
    int i, N=9;
    float kx[N], L=1;
    //TODO C Code for calculating frequencies
    //Frequencies 
    for (int i = 0; i< N / 2; ++i)
        kx[i] = i / L;

    kx[N / 2] = 0.00;

    for (int i = 0; i < ((N / 2) - 1); ++i)
        kx[i + 1 + N / 2] = -kx[N / 2 - i - 1];

    for (i = 0; i < N; i++)
        printf("%f\n", kx[i]);
    
}