#include "cochlea.h"
#include <cmath>
using namespace std;

void run(int n_osc, int n_t, float fs, double* stimulus, double* X_t)
{
    //SETUP
   
    setup(n_osc, n_t, fs, stimulus, X_t);

    //ITARATE
    for(int n = 0; n < n_t ; n++)
    {
        update();
    }

    delete_all();

}
