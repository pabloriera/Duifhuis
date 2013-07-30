#include <iostream>
#include "duifhuis.h"
#include "mex.h"


extern void _main();

//------------------------------------------------------------------

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{    
    int        n_osc, n_t;
    float      fs;
    double*    X_t;
    double*    stimulus;
    
    //Check for proper number of arguments 
    
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:Duifhius:nargin",
                "Requires four input arguments.");
    }
    
    //PARAMTERS READING
    n_osc = (int)*mxGetPr(prhs[0]);
    n_t = (int)*mxGetPr(prhs[1]);
    fs = (float)*mxGetPr(prhs[2]);
    stimulus = mxGetPr(prhs[3]);
    
    cout <<  "N_OSC = " << n_osc << endl;
    cout <<  "N_T = " << n_t << endl;          
  
    /*  create a C pointer to a copy of the output matrix */
    plhs[0] = mxCreateDoubleMatrix( (mwSize) n_osc*2 , (mwSize) n_t, mxREAL);
    X_t = mxGetPr(plhs[0]);
    
    run(n_osc,n_t,fs,stimulus,X_t);    
}
