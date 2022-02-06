// TODO: Using FFTW to denoise a noisy data

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>       // for generating random numbers (different each time)
#include <stdlib.h>     // rand(), srand()

int main()
{
    float dt = 0.001;
    int N=1001, i;
    double *f, f_clean[N], t[N];
    // FFTW DataTypes
    fftw_complex *f_hat;
    fftw_plan plan;

    // FFTW way of allocating array (There's no such issue in FORTRAN!!)
    f = (double *) fftw_malloc(sizeof(double) * N);
    // f_clean = (double *) fftw_malloc(sizeof(double) * N);
    f_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    srand(time(0));

    // the function
    for(i=0; i<N; i++)
    {
        t[i] = i*dt;
        f_clean[i] = sin(2*M_PI*50*t[i]) + sin(2*M_PI*150*t[i]) + sin(2*M_PI*200*t[i]);
        f[i] = f_clean[i] + rand();
    }


    //? Creating a file in C to export data to that file for plotting
    FILE *t_dat, *f_clean_dat, *f_noisy_dat, *f_hat_dat;
    t_dat = fopen("t.dat", "w");
    f_clean_dat = fopen("f_clean.dat", "w");
    f_noisy_dat = fopen("f_noisy.dat", "w");
    f_hat_dat = fopen("f_hat.dat", "w");

    for(i=0;i<N;i++)
    {
        fprintf(t_dat, "%f\n", t[i]);
        fprintf(f_clean_dat, "%f\n", f_clean[i]);
        fprintf(f_noisy_dat, "%f\n", f[i]);
    }




    // // FFTW Plan
    // plan = fftw_plan_dft_r2c_1d(N, f, f_hat, FFTW_ESTIMATE);

    // // Executing Plan
    // fftw_execute_dft_r2c(plan, f, f_hat);

    // Termination
    // fftw_destroy_plan(plan);
    fftw_free(f); fftw_free(f_hat);
    fclose(f_clean_dat); // files should be closed when done
    fclose(f_noisy_dat);
    fclose(f_hat_dat);
}