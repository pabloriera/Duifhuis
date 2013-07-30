#include <iostream>
#include <vector>
#include "duifhuis.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
    double* X_t ;
    int n_t = 8000;
    int n_osc = 200;
    float fs = 400000.0;

    double* stimulus;
    stimulus = new double[n_t+1];

    for(int n = 0; n < n_t; n++)
        stimulus[n] = 1*sin(6.28 * n / fs * 16000);

    mat aux(stimulus,n_t+1,1,false);
    aux.save("stimulus.txt",raw_ascii);

    X_t = new double[n_t*n_osc*2];

    run(n_osc, n_t, fs, stimulus,X_t);

    mat aux2(X_t, n_osc*2,n_t, false);
    aux2.save("X_t.txt", raw_ascii);

    delete[] X_t;

    return 0;
}
