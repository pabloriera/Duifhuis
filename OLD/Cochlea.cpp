#include <iostream>
#include "Cochlea.h"
#include <math.h>
#include "rk4.h"

using namespace std;

double SUM(double *array,int n);

Cochlea::Cochlea()
{
    
}

void Cochlea::setup(int _n)
{
    n = _n;
    n_eq = (n+1)*2;
    
    bm_length = cochleaLength - helicotremaWidth;
    bm_width  = scalaWidth;
    bm_mass   = bmMass * bmImpedanceFactor;
    
    dx = bm_length / n;
    
    for(int i = 0;i<n;i++)
        x[i] = dx * i;                
                
    Y  = new double[n+1];
    V  = new double[n+1];
    dY = new double[n+1];
    dV = new double[n+1];
    
    g  = new double[n+1];
    q  = new double[n+1];
    k  = new double[n+1];
    b  = new double[n+1];
    
    set_stiffness_damping();
    
    gaussian_elimination_init();
    
    X = new double[n_eq];
    for(int i = 0;i<n_eq;i++)
        X[i]=0;    
}

void Cochlea::gaussian_elimination_init()
{
    Scala_Area = scalaWidth * scalaHeight;
    Asq = 2 * rho * bm_width / (bm_mass * Scala_Area) * dx * dx;
    
    AME = 1 + (gam0 * Asq);
    ASC = 2 + Asq;
    AHT = 1 + Asq + sqrt(Asq);
    
    //frequency independent matching impedance when Qn=0,5
    phi = 2 * rho * dx / Scala_Area * sqrt(sdivm[n]);	// sdivm is already the stiffness divided by the mass, so the m is removed from the formula
    
    Asq0 = gam0 * Asq;
    
    if (useApexShortcut)
    {
        AHT = ASC;
        phi = 0;
    }
    
    /* 	 Calculation of coefficients necessary to solve the tridiagonal system.
     * (First the tridiagonal matrix is transformed into a lower triangle
     * matrix.
     * Furthermore we use the fact that matrix A is represented by only
     * three different diagonal elements (AME,ASC and AHT) and all sub-
     * and super diagonal elements are equal to -1.)
     */
    
    b[n] = 1 / AHT;
    for(int i = n-1;i>0;i--)
        b[i] = 1 / (ASC - b[i+1]);
    
    b[0] = 1 / (AME - b[1] );
    
}
void Cochlea::update()
{
    rk4(&equations, X, n_eq, t, dt);
    t += dt;
}

void Cochlea::equations(int n, double X[], double dX[], double t )
{    
    /*
     * n + 1    = number of oscilators
     * v[]  = vector with current state for in out
     * dv[] = derivative vector
     * t    = time
     */
    
    /* Systems equations
     * dY/dt = V;
     * dV/dt = q - g;
     * Y = v[even]
     * V = v[odd]
     */
    
//FORCE or STIMULUS
    
    double n_t = t / dt;
    int n_t_1  = ceil(n_t)+1;
    int n_t_2  = n_t_1 + 1;
    double ee = n_t - n_t_1;
    F0 = stimulus[n_t_1]*(1-ee) + stimulus[n_t_2]*ee;
    
//ASSIGN VARIABLE X
    
    Y = X;
    V = X + n + 1;
    
//COMPUTE g = d*V + s*Y
    
    //LINEAR CASE
    g[0] = d_m_factor * V[0] + s_m_factor * Y[0];
    
    for(int i=1;i<=n;i++)
    {
        g[i] = V[i] * ddivm[i] +Y[i] * sdivm[i];
    }
    
//COPMUTE q solving Aq = (k?);
    
    U_helicotrema = -stapesArea * V[0] - (SUM(V,n+1)-V[0]) * dx * bm_width;
    
    k[n] = -Asq * g[n] - phi * U_helicotrema;
    
    for(int i=n-1;i>0;i--)
        k[i] = -Asq * g[i] + k[i+1] * b[i+1];
    
    k[0] = -Asq0 * ( p0x * F0 + g[0] + r_Xtr0 * Y[0] ) +k[1] * b[1];
    q[0] = -k[0] * b[0];
    
    for(int i=1;i<=n;i++)
        q[i] = ( q[i-1] - k[i] ) * b[i];
    
    
//COMPUTE derivates of vector field
    
    dY[0] = V[0];
    dV[0] = m0_RK4 * (q[0] - g[0] - p0x * F0 - r_Xtr0 * Y[0]) * dt;
    
    for(int i=1;i<=n;i++)
    {
        dY[i] = V[i];
        dV[i] = q[i]-g[i];
        
    }
    
//RETURN VARIABLE dX
    
    for(int i=0;i<=n;i++)
    {
        dX[i] = dY[i];
        dX[i+n] = dV[i];
    }
    
    
    return;
}

double SUM(double *array,int n)
{
    double sum=0;
    
    for(int i=0;i<n;i++)
        sum += array[i];
    
    return sum;
}

void Cochlea::set_stiffness_damping()
{
    sdivm  = new double[n+1];
    ddivm  = new double[n+1];
    f_resonance  = new double[n+1];
    omega  = new double[n+1];
    
    if (useGreenwoodMap)
    {
        for(int i = 0;i<n;i++)
            f_resonance[i] = f_base_Greenwood_map * pow(10,(-kappa_Greenwood_map * x[i])) - apex_corr_Greenwood_map;
    }
    else
        for(int i = 0;i<n;i++)
            f_resonance[i] = f_base_exp_map * pow(10, ( -kappa_exp_map * x[i] ));
    
    for(int i = 0;i<n;i++)
        omega[i] = 2 * PI * f_resonance[i];
    
    if (useConstantQ)
    {
        QualityFactor = ConstantQ;
        for(int i = 0;i<n;i++)
            ddivm[i] = omega[i] / QualityFactor;
    }
    else
    {
        /*
         * Use a square root profile for Q
         * Make Q equal to 0.5 at the apex
         * in order to get an optimal, frequency independent,
         * matching impedance at that point.
         * Q(i) = 0.5 + const * SQRT( f(i) - f(n) )
         * const is set so Q(i)>20 for f(i) > 4 kHz
         * original version JK		QualityFactor = 0.50 + 0.30 * DSQRT(f_resonance - f_resonance(n))
         * original HD       		QualityFactor = 0.50*DSQRT(f_resonance/f_resonance(n))
         *		modified HD 21 aug 2005; do not use Q-factor as primary parameter
         *  use omega and sdivm, then Q and ddivm
         *! test 2009: reduce ddivm by factor 2
         *		ddivm=1.100*omega**(1.0/3.0) original
         */
        for(int i = 0;i<n;i++)
            ddivm[i]=100*pow(omega[i],(1.0/3.0)) / 2.0;
    }
    
    /* Definition of Quality factor ("Classical Dynamics of particles and systems" by M. Thornton, page 128):
     * Q = omega * mass / damping
     * ---> ddivm in program already scaled by mass --->
     * ddivm = omega / QualityFactor*/
    for(int i = 0;i<n;i++)
        sdivm[i] = omega[i]*omega[i];
}

void Cochlea::parameters()
{
    
    cochleaLength                  = 36e-3;
    bmMass                         = 0.5;
    bmImpedanceFactor              = 1;
    scalaWidth                     = 1e-3;
    scalaHeight                    = 1e-3;
    helicotremaWidth               = 1e-3;
    rho                            = 1000;
    Q_constant_default             = 20;
    f_base_exp_map                 = 22507;
    kappa_exp_map                  = 65.1;
    f_base_Greenwood_map           = 20682;
    kappa_Greenwood_map            = 60;
    apex_corr_Greenwood_map        = 140.59;
    stapesArea                     = 3e-6;
    EardrumArea                    = 60e-6;
    MiddleEarResonanceFrequency    = 2e3;
    MiddleEarQualityFactor         = 0.4;
    SpecificAcousticImpedanceOfAir = 415;
    middleEarTransformer           = 30;
    damping_coupler                = 140e5;
    mass_coupler                   = 43.4e2;
    stiffness_coupler              = 1/2.28e-11;
    p0                             = 2e-5;
    
    /*
     * cochleaLength     = [36d-3, 18.5d-3]  !m      **HUMAN set to 35 + 1mm for apex [Dec09]
     * bmMass(2)             = [0.5d0, 0.25d0]	!kg/m2
     * bmImpedanceFactor(2)	= [1d0, 0.75d0]     !dimensionless
     * !bmInpedanceFactor: scaling factor which enables overall parameter adjustments. Default: := 1
     * scalaWidth(2)         = [1d-3, 0.5d-3]	!m
     * scalaHeight(2)        = [1d-3, 0.5d-3]	!m
     * helicotremaWidth(2)   = [1d-3, 0.5d-3]	!m  note: also used for helicotremalength
     * rho                   = 1d3		        !kg/m3
     * rho: average density of cochlear fluid, approx. equal to density of water
     * Q_constant_default    = 20d0              !default Q-value for cochlear tuning curve, DIM 0
     * bm_mass	    ! kg/m2
     * bm_width	! m
     * bm_length	! m
     *
     * Frequency map parameters, x measured from stapes
     * Exponential Map:  f(x) = f_base * 10 ^ (-kappa x)
     * f_base_exp_map (2)    = [22507d0, 38346d0]	!Hz
     * kappa_exp_map (2)     = [65.1d0, 113.5d0]	    !1/m
     *
     * Greenwood Map:    f(x) = f_base * 10 ^ (-kappa x) - f_apex  *****was formulation before 2010
     **************************************************************************************************Modified
     * The original formula is f = 165.4 ( 10 ^ (kappa x1) - k ) with some discussion about k, but kappa = 60 /m
     *   where x1 is measured from apex. Using x = 35 - x1, gives:
     *   f = 165.4 ( 10 ^ (60*(0.035-x))- k ) = 165.4 * ( 10 ^2.1 * 10 ^ (-kappa x) - k)
     *     = 165.4 * 10^2.1 ( 10 ^ (-kappa x) - k * 10^(-2.1) )
     *     = 165.4 * 10^2.1 ( 10 ^ (-kappa x)) - k * 165.4
     *   The apex value = 165.4 ( 10^0 - k) = 165.4 (1 - k) !!! => zero for k=1; 24.81 for k=0.85; and 33.08 for k=0.8 .
     *   value used originally by us is: k = 0.85
     *   Likewise, f_base= 165.4 (10^2.1 - k), => 20674 for k=0.9; 20682 for k=0.85; 20690 for k=0.8
     *   This has been the basis for the value: f_base_Greenwood_map.
     *   The value originally called f_apex was just k x 165.4! Formula is in principle OK, but term f_apex should be called apex_corr
     *   and values have to be adjusted to the correct values.  Effect correction appears marginal
     *
     * f_base_Greenwood_map(2)     = [20682d0,  43765d0]	!Hz
     * kappa_Greenwood_map (2)     = [60d0,     113.5d0]	!1/m (=2.1/bmlength)
     * apex_corr_Greenwood_map (2) = [140.59d0, 297.5d0]	!Hz not the proper term. It is 165.4 x k-term
     *
     * Middle Ear parameters
     * stapesArea(2)                 = [3d-6, 0.81d-6]   !m2
     * EardrumArea(2)                = [60d-6, 23.9d-6]  !m2
     * MiddleEarResonanceFrequency(2)= [2d3, 4d3]        !Hz
     * MiddleEarQualityFactor        = 0.4d0             !dimensionless
     * SpecificAcousticImpedanceOfAir = 415d0            !Ns/m3        = Pa/(m/s)
     * middleEarTransformer(2)       = [30d0, 39.3d0]    !dimensionless
     * damping_coupler(2)            = [140d5, 280d5]	!Ns/m5 =kg/sm4= Pa/(m3/s)
     * mass_coupler(2)               = [43.4d2, 21.7d2]	!kg/m4
     * stiffness_coupler             = 1d0/2.28d-11      !kg/s2m4
     *
     * p0    = 2d-5  !Pa = N/m2
     */
       
}

