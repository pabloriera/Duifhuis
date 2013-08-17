#include <cmath>

//#define M_PI 3.14159265358979

class Cochlea_t
{public:
    Cochlea_t(){}
    
    void gaussian_elimination_init();
    void equations(int n, double X[], double dX[], double t );
    void set_arrays();
    void constants();
    void parameters();
    void rk4();
    
    void setup(int _n_osc,int _n_t, float fs, double* _stimulus, double* Y_t, double* V_t);
    void update();
    void delete_all();
        
    double *F1, *F2, *F3, *F4, *Xtemp;
    
    double* Y;
    double* V;
    double* dY;
    double* dV;
    
    double* g;
    double* k;
    double* q;
    double* b;
    
    double* sdivm;
    double* ddivm;
    
    double* f_resonance;
    double* omega;
    
    double* X;
    double* x;
    
    double* stimulus;
    
    double* Y_t;
    double* V_t;
    
    double t, dt;
    
    int N, n_eq;
    int c,n_t;
    
    double AME,ASC,AHT,Asq;
    double Scala_Area, scalaWidth, scalaHeight;
    double phi, Asq0, gam0;
    double rho,bm_mass,bm_width;
    double dx, cochleaLength, bm_length;
    double bmMass, bmImpedanceFactor;
    double helicotremaWidth;
    double F0;
    double dv1,dv3,dvy2,dy3;
    
    bool useGreenwoodMap;
    bool useConstantQ;
    bool useApexShortcut;
    bool useEarCanalCoupler;
    
    double QualityFactor,ConstantQ,Q_constant_default;
    double p0,p0x;
    double r_Xtr0,m0_RK4;
    
    double f_base_exp_map;
    double kappa_exp_map;
    double f_base_Greenwood_map;
    double kappa_Greenwood_map;
    double apex_corr_Greenwood_map;
    double stapesArea;
    double EardrumArea;
    double MiddleEarResonanceFrequency;
    double MiddleEarQualityFactor;
    double SpecificAcousticImpedanceOfAir;
    double middleEarTransformer;
    double damping_coupler;
    double mass_coupler;
    double stiffness_coupler;
    
    double g0_factor;
    double q0_factor;
    double Y0_factor;
    double d_m_factor;
    double s_m_factor;
    
    double U_helicotrema;
    
};

