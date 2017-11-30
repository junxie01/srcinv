
#ifndef mc_hpp
#define mc_hpp

#include <stdio.h>
#include "aux.hpp"

#endif /* mc_hpp */


/********************************************************************************//**
* Base class for inversion.
************************************************************************************/
class Inversion
{
public:

    Input in;
};

/********************************************************************************//**
* Derived class for Monte Carlo inversion methods.
************************************************************************************/
class MonteCarlo: public Inversion
{
public:
    
    /* Model vectors, current and testing. */
    double *q;              /**< Current model vector. */
    double *q_new;          /**< Test model vector. */
    
    /* Writing samples to a file. */
    void write_sample(FILE *pfile,                  /**< Pointer to open file. */
                      double misfit,                /**< Misfit value. */
                      int iteration,                /**< Current iteration. */
                      int multiplicity              /**< Multiplicity of the sample. */
    );
};

/********************************************************************************//**
* Derived class for Metropolis-Hastings inversion.
************************************************************************************/
class MH: public MonteCarlo
{
public:
    /* Proposal and misfit functions. ---------------------------------------------*/
    void propose();             /**< Proposal based on Hamiltonian mechanics. */
    double likelihood();        /** Likelihood function. */
    double misfit();            /** Posterior misfit. */
    
    /* Constructor and destructor. ------------------------------------------------*/
    MH(const char *input_filename);
    ~MH();
};


/********************************************************************************//**
* Derived class for Hamiltonian Monte Carlo inversion.
************************************************************************************/
class HMC: public MonteCarlo
{
public:
    /* HMC-specific fields. -------------------------------------------------------*/
    double *p;                  /**< Current momentum vector. */
    double *p_new;              /**< Test momentum vector. */
    double *m;                  /**< Mass matrix diagonal elements. */
    
    /* Right-hand sides of Hamilton equations. ------------------------------------*/
    void dHdp(                  /** dH/dp. Needs to be adjusted as needed. */
              double *p,        /**< Momentum vector. */
              double *rhs       /**< Output right-hand side, i.e. dH/dp. Must be allocated. */
    );
    
    void dHdq(                  /** dH/dq. Needs to be adjusted as needed. */
              double *q,        /**< Model (position) vector. */
              double *rhs       /**< Output right-hand side, i.e. dH/dq. Must be allocated. */
    );
    
    /* Leap-frog integration of Hamilton's equations with initial positions and momenta. */
    void leap_frog();
    
    /* Proposal and energy functions. ---------------------------------------------*/
    void propose();             /**< Proposal based on Hamiltonian mechanics. */
    double energy();            /** Total energy for Hamiltonian Monte Carlo. */
    double potential_energy();  /** Potential energy, i.e. posterior. */

    /* Constructor and destructor. ------------------------------------------------*/
    HMC(const char *input_filename);
    ~HMC();
};

