#include <iostream>
#include <cmath>
#include <fstream>
#include <fftw3.h>

using namespace std;
const double Pi = 3.141592653589;
const int N = 64;
const double L = 5.0;

void Derivative(double& num1, double& num2, double J) {
    //takes two derivatives
    num1 = -1 * pow(J, 2)*num1;
    num2 = -1 * pow(J, 2)*num2;
}

double Function(double X) { // Initial function
    double value = 1.0 / cosh(3 * (X - M_PI));//sin(X) + sin(3*X);
    return value;
}

void FFT(double *in, double *Result) {

    double *input, *outI, *kx;
    input = new double[N];
    outI = new double[N];
    kx = new double[N];

    for (int i = 0; i<N; i++) {
        input[i] = in[i];
    }

    //Declaring transformed matricies;
    fftw_complex *out;
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N));

    //Plans for execution
    fftw_plan p;

    p = fftw_plan_dft_r2c_1d(N, input, out, FFTW_ESTIMATE);

    fftw_execute(p);

    //----This part calculate frequencies--------//
    //Frequencies 
    for (int i = 0; i< N / 2; ++i) {
        kx[i] = i / L;
    }

    kx[N / 2] = 0.00;

    for (int i = 0; i < ((N / 2) - 1); ++i) {
        kx[i + 1 + N / 2] = -kx[N / 2 - i - 1];
    }

    //Derivative
    for (int i = 0; i<N; i++) {
        //the frequencies is passed to this function to calculate
        //the derivatives
        Derivative(out[i][0], out[i][1], kx[i]);
    }

    fftw_plan pI;
    //Execution of Inverse FFT

    pI = fftw_plan_dft_c2r_1d(N, out, outI, FFTW_PRESERVE_INPUT);

    fftw_execute(pI);

    for (int i = 0; i<N; i++) {//Dividing by the size of the array
        Result[i] = (1.0 / N)*outI[i];
    }

    fftw_free(outI); fftw_free(out);
    fftw_free(input);
    fftw_destroy_plan(pI);  fftw_destroy_plan(p);
}

int main()
{
    //Creating and delcaring function variables
    double *initial, *w1, *w2, *w3, *w4, *y, *holder1, *holder2, *holder3;
    double dx = 2 * Pi / N, x = 0, t = 0, dt = 0.01;
    double h = pow(dx, 2), tmax = 100 * h;

    //creating arrays for the RK4 procedure
    initial = new double[N];
    w1 = new double[N]; w2 = new double[N]; y = new double[N];
    w3 = new double[N]; w4 = new double[N]; holder1 = new double[N];
    holder2 = new double[N]; holder3 = new double[N];
    double *resultArray = new double[N];
    double *resultnext = new double[N];
    int j = 1;

    // The stability parameter are declared but not used
    float alpha = 0.1, CFL = alpha * dt / (10 * pow(dx, 2));

    //Output files
    ofstream sendtofileINITIAL("HeatdataINITIAL.dat");
    ofstream sendtofile5("Heatdata5.dat");
    ofstream sendtofile25("Heatdata25.dat");
    ofstream sendtofile85("Heatdata85.dat");
    ofstream sendtofileFINAL("HeatdataFINAL.dat");

    //Initial data
    for (int i = 0; i < N; i++) {
        y[i] = Function(x);
        sendtofileINITIAL << i * 2 * Pi / N << " " << y[i] << endl;
        x += dx;
    }

    //RK4
    // y[i+1] = y[i] + (h/2)*(w1 + 2*w2 + 2*w3 + w4)
    // where the w1,w2,w3,w4 is the standard RK4 calculations
    // w1= f(y[i]), w2= f(y[i]+(h/2)*w1[i])
    // w3= f(y[i] + (h/2)*w2[i]) and w4 = f(y[i]+h*w3[i])
    while (t <= tmax) {

        FFT(y, resultArray);

        for (int i = 0; i < N; i++) {    //Calculating w1
            w1[i] = resultArray[i];
        }

        for (int i = 0; i < N; i++) {
            holder1[i] = y[i] + (h / 2)*w1[i];
        }

        FFT(holder1, resultArray);

        for (int i = 0; i < N; i++) {   //Calculating w2
            w2[i] = resultArray[i];
        }

        for (int i = 0; i < N; i++) {
            holder2[i] = y[i] + (h / 2)*w2[i];
        }

        FFT(holder2, resultArray);

        for (int i = 0; i < N; i++) {   //Calculating w3
            w3[i] = resultArray[i];
        }

        for (int i = 0; i < N; i++) {
            holder3[i] = y[i] + h * w3[i];
        }

        FFT(holder3, resultArray);

        for (int i = 0; i < N; i++) {   //Calculating w4
            w4[i] = resultArray[i];
        }

        for (int i = 0; i < N; i++) {
            y[i] += (h / 6.0)*(w1[i] + 2 * w2[i] + 2 * w3[i] + w4[i]);

            //Outputing data at certain time periods
            if (j == 5)
                sendtofile5 << i * 2 * Pi / N << " " << y[i] << endl;

            if (j == 25)
                sendtofile25 << i * 2 * Pi / N << " " << y[i] << endl;

            if (j == 85)
                sendtofile85 << i * 2 * Pi / N << " " << y[i] << endl;

            if (j == 100)
                sendtofileFINAL << i * 2 * Pi / N << " " << y[i] << endl;
        }
        j++;// counter
        t += h;// time counter
    }
    sendtofileINITIAL.close();
    sendtofile5.close();
    sendtofile25.close();
    sendtofile85.close();
    sendtofileFINAL.close();
    return 0;
}