/* Runge Kutta integrator from numerical recipies plus improvements */
/* void *deri(int n,double h[],double D[],double t);  */
/* function argument not tested yet */
#include "rk4.h"

void rk4(void deri(int , double* [], double* [], double ), double h[], int n, double t, double dt)
{

    int i;
    double dt2, dt6;

    double *k1,*k2,*k3,*k4,*h0;

    k1 = new double[n];
    k2 = new double[n];
    k3 = new double[n];
    k4 = new double[n];
    h0 = new double[n];


    dt2=dt/2.0;
    dt6=dt/6.0;
/*
    cout << "h_in ";
    for (i = 0 ; i<n; i++)
        cout << h[i] << " ";
    cout << endl;*/

    for (i = 0 ; i<n; i++)
        h0[i] = h[i];

    deri(n,&h0,&k1,t);
    for (i =0 ; i<n ; i++)
        h0[i] = h[i]+dt2*k1[i];

    deri(n,&h0,&k2,t+dt2);
    for (i =0 ; i<n ; i++)
        h0[i] = h[i]+dt2*k2[i];

    deri(n,&h0,&k3,t+dt2);
    for (i =0 ; i<n ; i++)
        h0[i] = h[i]+dt*k3[i];

    deri(n,&h0,&k4,t+dt);
    for (i = 0 ; i<n ; i++)
        h0[i] = h[i]+dt*k4[i];


    for (i =0; i<n ; i++)
            h[i] = h[i] + dt6*(2.*(k2[i]+k3[i])+k1[i]+k4[i]);

/*    cout << "h_out ";
    for (i = 0 ; i<n; i++)
        cout << h[i] << " ";
    cout << endl;*/

    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] h0;

}
