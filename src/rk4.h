#include <iostream>
#include <vector>

using namespace std;

void rk4(double X[], int nX, double t, double dt, void (*derivsRK)(int nX, double X[], double deriv[], double t ) );
