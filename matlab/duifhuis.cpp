#include <iostream>
#include "mex.h"
#include "src/cochlea.h"

using namespace std;

extern void _main();

//------------------------------------------------------------------

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int        n_osc, n_t;
    float      fs, rho;
    double dv1,dv3,dvy2,dy3;
    double     *Y_t, *V_t;
    double*    stimulus;
    Cochlea_t  cochlea;
    
    //Check for proper number of arguments
    
    if (nrhs < 2 || nrhs > 5) {
        mexErrMsgIdAndTxt("MATLAB:Duifhius:nargin",
                "Requires minimun two input arguments. Stimulus, fs, n_t, n_osc,rho, damping_coeffs" );
    }
    
    //Input parameters
    if (nrhs >= 2)
    {
        stimulus = mxGetPr(prhs[0]);
        n_t = mxGetM(prhs[0]);
        if(n_t==1)
            n_t = mxGetN(prhs[0]);
        //cout << "n_t = " << n_t << endl;
        
        fs = (float)*mxGetPr(prhs[1]);
    }
    
    if (nrhs >= 3)
        n_osc = (int)*mxGetPr(prhs[2]);
    else
    {
        cout << "n_osc = 300" << endl;
        n_osc = 300;
    }
    
    if (nrhs >= 4)
        rho = (float)*mxGetPr(prhs[3]);
    else
    {
        cout << "rho = 1000" << endl;
        rho = 1000;
    }
    
    if (nrhs >= 5)
    {
        int nfields = mxGetNumberOfFields(prhs[4]);
        
        if(nfields == 4)
        {
            double *ptr = mxGetPr(mxGetField(prhs[4], 0, "dv1"));
            dv1 = *ptr;
            ptr = mxGetPr(mxGetField(prhs[4], 0, "dv3"));
            dv3 = *ptr;
            ptr = mxGetPr(mxGetField(prhs[4], 0, "dvy2"));
            dvy2 = *ptr;
            ptr = mxGetPr(mxGetField(prhs[4], 0, "dy3"));
            dy3 = *ptr;
        }
        else
            mexErrMsgIdAndTxt("MATLAB:Duifhius:vargin",
                    "Requires a structure with these fields, dv1, dv3, dy1,dy3" );
    }
    else
    {
        dv1 = 1;
        dv3 = 0;
        dvy2 = 0;
        dy3 = 0;
    }
    
    
    /*  create a C pointer to a copy of the output matrix */
    plhs[0] = mxCreateDoubleMatrix( (mwSize) n_osc , (mwSize) n_t, mxREAL);
    Y_t = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix( (mwSize) n_osc , (mwSize) n_t, mxREAL);
    V_t = mxGetPr(plhs[1]);
    
    
    //Start cochlea
    cochlea.setup(n_osc, n_t, fs, stimulus, Y_t, V_t);
    
    //Pass extra input parameters
    cochlea.rho = rho;
    cochlea.dv1 = dv1;
    cochlea.dv3 = dv3;
    cochlea.dvy2 = dvy2;
    cochlea.dy3 = dy3;   
    
    
    cochlea.gaussian_elimination_init();
    
    //Iterate
    for(int n = 0; n < n_t ; n++)
    {
        cochlea.update();
    }
}