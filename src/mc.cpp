#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mc.hpp"
#include "aux.hpp"


/********************************************************************************//**
* Derived class for Monte Carlo inversion methods.
************************************************************************************/

/* Write a sample to an open file. ------------------------------------------------*/
void MonteCarlo::write_sample(FILE *pfile, double misfit, int iteration, int multiplicity)
{
    if (iteration==0) fprintf(pfile,"%d %d\n",in.dim,in.n_samples+1);
    
    for (int i=0; i<in.dim; i++) fprintf(pfile,"%lg ",in.scale_q[i]*q[i]);
    fprintf(pfile,"%lg %d\n",misfit,multiplicity);
}

/********************************************************************************//**
* Derived class for Metropolis-Hastings inversion.
************************************************************************************/

/* Constructor. -------------------------------------------------------------------*/
MH::MH(const char *input_filename)
{
    /* Read input. ----------------------------------------------------------------*/
    in.read(input_filename);
    
    /* Allocate memory. -----------------------------------------------------------*/
    q=new double[in.dim];
    q_new=new double[in.dim];

    /* Initialise random number generator. ----------------------------------------*/
    srand (time(0));
    
    /* Initialise models ----------------------------------------------------------*/
    for (int i=0; i<in.dim; i++)
    {
        q_new[i]=randn(in.mean_q[i],in.sigma_q[i]); /* Radom start around prior mean. */
        q[i]=q_new[i];
    }
}

/* Destructor. -------------------------------------------------------------------*/
MH::~MH()
{
    if (q) delete[] q;
    if (q_new) delete[] q_new;
}

/* Propose test model based on prior. ---------------------------------------------*/
void MH::propose()
{
    for (int i=0; i<in.dim; i++) q_new[i]=randn(in.mean_q[i],in.sigma_q[i]);
}

/* Likelihood function for Metropolis Hastings. -----------------------------------*/
double MH::likelihood()
{
    /* Local variables. */
    double chi=0.5*in.C;
    
    /* Compute likelihood. */
    for (int i=0; i<in.dim; i++)
    {
        chi+=(q_new[i]-in.mean_q[i])*in.B[i];
        for (int j=0; j<in.dim; j++) chi+=0.5*(q_new[i]-in.mean_q[i])*(q_new[j]-in.mean_q[j])*in.A[i][j];
    }
    return chi;
}

/* Posterior misfit for Metropolis Hastings. ---------------------------------------*/
double MH::misfit()
{
    /* Local variables. */
    double chi=0.5*in.C;
    
    /* Compute likelihood. */
    for (int i=0; i<in.dim; i++)
    {
        chi+=(q_new[i]-in.mean_q[i])*in.B[i];
        for (int j=0; j<in.dim; j++) chi+=0.5*(q_new[i]-in.mean_q[i])*(q_new[j]-in.mean_q[j])*in.A[i][j];
    }
    
    /* Add prior. */
    for (int i=0; i<in.dim; i++) chi+=0.5*(q_new[i]-in.mean_q[i])*(q_new[i]-in.mean_q[i])/(in.sigma_q[i]*in.sigma_q[i]);
    
    return chi;
}



/********************************************************************************//**
* Derived class for Hamiltonian Monte Carlo inversion.
************************************************************************************/

/* Constructor. -------------------------------------------------------------------*/
HMC::HMC(const char *input_filename)
{
    /* Read input. ----------------------------------------------------------------*/
    in.read(input_filename);
    
    /* Allocate memory. -----------------------------------------------------------*/
    p=new double[in.dim];
    p_new=new double[in.dim];
    q=new double[in.dim];
    q_new=new double[in.dim];
    m=new double[in.dim];
    
    /* Initialise random number generator. ----------------------------------------*/
    srand (time(0));
    
    /* Add the prior. -------------------------------------------------------------*/
    for (int i=0; i<in.dim; i++) in.A[i][i]+=1.0/(in.sigma_q[i]*in.sigma_q[i]);
    
    /* Initialise mass matrix. ----------------------------------------------------*/
    for (int i=0; i<in.dim; i++) m[i]=in.A[i][i]+in.hmc_reg;
    
    /* Initialise models and momenta. ---------------------------------------------*/
    for (int i=0; i<in.dim; i++)
    {
        p_new[i]=randn(0.0,sqrt(m[i]));             /* Momentum according to mass matrix. */
        q_new[i]=randn(in.mean_q[i],in.sigma_q[i]); /* Radom start around prior mean. */
        q[i]=q_new[i];
    }
}

/* Destructor. -------------------------------------------------------------------*/
HMC::~HMC()
{
    if (p) delete[] p;
    if (p_new) delete[] p_new;
    if (q) delete[] q;
    if (q_new) delete[] q_new;
    if (m) delete[] m;
}

/* Right-hand sides of Hamilton equations. ----------------------------------------*/
void HMC::dHdp(double *p_in, double *rhs)
{
    for (int i=0; i<in.dim; i++) rhs[i]=p_in[i]/m[i];
}

void HMC::dHdq(double *q_in, double *rhs)
{
    /* March through components. */
    for (int k=0; k<in.dim; k++)
    {
        rhs[k]=in.B[k];
        for (int i=0; i<in.dim; i++) rhs[k]+=(q_in[i]-in.mean_q[i])*in.A[i][k];
    }
}

/* Leap-frog integration of Hamilton's equations. ---------------------------------*/
void HMC::leap_frog()
{
    /* Local variables and setup. -------------------------------------------------*/
    
    int it, nt;
    FILE *pfile;
    double* out=new double[in.dim];
    double* p_half=new double[in.dim];
    
    if (in.hmc_output_trajectory) pfile=fopen("../output/trajectory.txt","w");
    
    /* Determine length of integration. -------------------------------------------*/
    
    if (!strcmp(in.hmc_length,"random"))
    {
        nt=(int)randf(0.0,in.hmc_nt);
    }
    else
    {
        nt=in.hmc_nt;
    }
    
    /* March forward. -------------------------------------------------------------*/
    
    dHdq(q,out);
    
    for (it=0; it<nt; it++)
    {
        /* Some output. */
        if (in.hmc_output_trajectory)
        {
            for (int i=0; i<in.dim; i++) fprintf(pfile,"%lg ",q[i]);
            for (int i=0; i<in.dim; i++) fprintf(pfile,"%lg ",p[i]);
            fprintf(pfile,"\n");
        }
        
        /* First half step in momentum. */
        for (int i=0; i<in.dim; i++) p_half[i]=p[i]-0.5*in.hmc_dt*out[i];
        
        /* Full step in position. */
        dHdp(p_half,out);
        for (int i=0; i<in.dim; i++) q_new[i]=q[i]+in.hmc_dt*out[i];
        
        /* Second half step in momentum. */
        dHdq(q_new,out);
        for (int i=0; i<in.dim; i++) p_new[i]=p_half[i]-0.5*in.hmc_dt*out[i];
       
        /* Update position and momentum. */
        for (int i=0; i<in.dim; i++)
        {
            p[i]=p_new[i];
            q[i]=q_new[i];
        }
    }

    if (in.verbose) printf("integration steps: %d\n",it);
    
    /* Clean up. ------------------------------------------------------------------*/
    
    if (in.hmc_output_trajectory) fclose(pfile);
    
    delete[] p_half;
    delete[] out;
}

/* Proposal based on the solution of Hamilton's equation. -------------------------*/
void HMC::propose()
{
    /* Draw random prior momenta. */
    for (int i=0; i<in.dim; i++) p[i]=randn(0.0,sqrt(m[i]));
}

/* Misfit for Hamiltonian Monte Carlo. --------------------------------------------*/
double HMC::energy()
{
    /* Local variables. */
    double H=0.5*in.C;
    
    for (int i=0; i<in.dim; i++)
    {
        H+=(q[i]-in.mean_q[i])*in.B[i];
        for (int j=0; j<in.dim; j++) H+=0.5*(q[i]-in.mean_q[i])*(q[j]-in.mean_q[j])*in.A[i][j];
    }
    
    /* Compute energy - momentum part. */
    for (int i=0; i<in.dim; i++) H+=0.5*p[i]*p[i]/m[i];
    
    return H;
}

/* Potential energy, i.e. posterior. ---------------------------------------------*/
double HMC::potential_energy()
{
    /* Local variables. */
    double U;
    
    /* Compute energy - model part. */
    U=0.5*in.C;
    
    for (int i=0; i<in.dim; i++)
    {
        U+=(q[i]-in.mean_q[i])*in.B[i];
        for (int j=0; j<in.dim; j++) U+=0.5*(q[i]-in.mean_q[i])*(q[j]-in.mean_q[j])*in.A[i][j];
    }
    
    return U;
}



