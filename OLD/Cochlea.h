#define PI 3.14159265358979

class Cochlea {
    
public:
    
    Cochlea();
    
    void setup(int _n);
    void update();
    void set_stiffness_damping();
    void parameters();
    void gaussian_elimination_init();
    void equations(int n, double X[], double dX[], double t );
    //void rk4(void deri(int , double [], double [], double ), double h[], int n, double t, double dt);
    
    ~Cochlea() { };
            
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
    
    double t, dt;
    
    int n, n_eq;
    
    double AME,ASC,AHT,Asq;
    double Scala_Area, scalaWidth, scalaHeight;
    double phi, Asq0, gam0;
    double rho,bm_mass,bm_width;
    double dx, cochleaLength, bm_length;
    double bmMass, bmImpedanceFactor;
    double helicotremaWidth;
    double F0;
    
    bool useGreenwoodMap;
    bool useConstantQ;
    bool useApexShortcut;
    
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
    
    double d_m_factor;
    double s_m_factor;
    double U_helicotrema;
    
};

