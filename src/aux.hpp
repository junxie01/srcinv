
#ifndef aux_hpp
#define aux_hpp

/*=================================================================================*/
/* Input. -------------------------------------------------------------------------*/
/*=================================================================================*/

class Input
{
public:
    
    /** General setup. */
    char method[10];                        /**< Inversion method (HMC, MH, LIN). */
    int dim;                                /**< Model space dimension. */
    bool verbose;                           /**< Write basic screen output. */
    
    /** Model parameters. */
    double *mean_q;                         /**< Prior means. */
    double *sigma_q;                        /**< Prior standard deviations. */
    double *scale_q;                        /**< Scalings of the model parameters. */
    
    // For source inversion: q[0-5]: moment tensor components, Mxx, Mxy, Mxz, Myy, Myz, Mzz;
    // q[6-8]: source coordinates, lon, lat, z; q[9]: origin time.
    
    /* Precomputed matrices for fast misfit & derivative evaluation. */
    double **A;
    double *B;
    double C;
    
    /** Monte Carlo-specific parameters. */
    int n_samples;                          /**< Maximum number of samples. */
    
    /**< HMC-specific parameters. */
    int hmc_nt;                             /**< HMC maximum number of trajectory time steps. */
    double hmc_dt;                          /**< HMC trajectory time increment. */
    double hmc_gamma;                       /**< HMC gravitational constant. */
    double hmc_reg;                         /**< HMC mass matrix regularisation. */
    bool hmc_output_trajectory;                 /**< Write trajectory to a file. */
    
    /** Read input from file. */
    void read(const char *input_filename);       /**< Read input from "input_file". */
    
    /** Constructor and destructor. */
    Input(const char *input_filename);      /**< Read input during construction. */
    Input();                                /**< Do-nothing constructor. */
    ~Input();
};

/*=================================================================================*/
/* Modelling parameters. ----------------------------------------------------------*/
/*=================================================================================*/

class parameters
{
public:
    
    /** Modelling parameters. */
    
    double q[20], sigma_q[20], mean_q[20];          /**< Inversion parameters (values, prior stdev, prior mean).
                                                     q[0-5]: moment tensor components, Mxx, Mxy, Mxz, Myy, Myz, Mzz,
                                                     q[6-8]: source coordinates, lon, lat, z,
                                                     q[9]: origin time,
                                                     q[10-19]: source-time function parameters. */
    
    double scale[20];                               /**< Scaling values for the above parameters. */
    
    /** Constructor and destructor. */
    
    parameters();         /**< Constructor. */
    ~parameters();        /**< Destructor. */
    parameters(const parameters &a);                            // Copy constructor.
    parameters &operator=(const parameters &a);                 // Assignment operator.

    /** Member functions. */
    
    void read_input(void);      /**< Fill values by reading INPUT/parameters_source.txt and INPUT/parameters_scaling.txt. */
};

/* Operators. ---------------------------------------------------------------------*/

parameters operator+(const parameters &a, const parameters &b);
parameters operator-(const parameters &a, const parameters &b);


/*=================================================================================*/
/* Little helpers. ----------------------------------------------------------------*/
/*=================================================================================*/

/** Double-valued, uniformly distributed random numbers. */
double randf(
             double min,    /**< Minimum value. */
             double max     /**< Maximum value. */
);

/** Double-valued, normally distributed random numbers. */
void randn(
           double mean,        /**< Mean. */
           double stdv,        /**< Standard deviation. */
           double *x1,         /**< Pointer to first random number. */
           double *x2          /**< Pointer to second random number. */
);

double randn(
             double mean,       /**< Mean. */
             double stdv        /**< Standard deviation. */
);

#endif /* aux_hpp */
