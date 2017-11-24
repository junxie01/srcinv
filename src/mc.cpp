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
    
    /* Regularisation / prior. ----------------------------------------------------*/
    for (int i=0; i<in.dim; i++) in.A[i][i]+=in.A[i][i]+0.5/(in.sigma_q[i]*in.sigma_q[i]);
    
    /* Initialise mass matrix. ----------------------------------------------------*/
    for (int i=0; i<in.dim; i++) m[i]=in.hmc_gamma*(in.A[i][i]+in.hmc_reg);
    
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
void HMC::leap_frog(bool output_trajectory)
{
    /* Local variables and setup. -------------------------------------------------*/
    
    int it;
    FILE *pfile;
    double angle1, angle2;
    double* out=new double[in.dim];
    double* p_half=new double[in.dim];
    double* p_init=new double[in.dim];
    double* q_init=new double[in.dim];
    
    if (output_trajectory) pfile=fopen("../output/trajectory.txt","a");
    
    /* Set initial values (current state). ----------------------------------------*/
    
    for (int i=0; i<in.dim; i++)
    {
        p_init[i]=p[i];
        q_init[i]=q[i];
    }
    
    /* March forward. -------------------------------------------------------------*/
    
    if (output_trajectory)
    {
        for (int i=0; i<2*in.dim; i++) fprintf(pfile,"0.0 ");
        fprintf(pfile,"\n");
    }
    
    dHdq(q_init,out);
    
    for (it=0; it<in.hmc_nt; it++)
    {
        /* Some output. */
        if (output_trajectory)
        {
            for (int i=0; i<in.dim; i++) fprintf(pfile,"%lg ",q_init[i]);
            for (int i=0; i<in.dim; i++) fprintf(pfile,"%lg ",p_init[i]);
            fprintf(pfile,"\n");
        }
        
        /* First half step in momentum. */
        for (int i=0; i<in.dim; i++) p_half[i]=p_init[i]-0.5*in.hmc_dt*out[i];
        
        /* Full step in position. */
        dHdp(p_half,out);
        for (int i=0; i<in.dim; i++) q_new[i]=q_init[i]+in.hmc_dt*out[i];
       
        /* Second half step in momentum. */
        dHdq(q_new,out);
        for (int i=0; i<in.dim; i++) p_new[i]=p_half[i]-0.5*in.hmc_dt*out[i];
       
        /* Update position and momentum. */
        for (int i=0; i<in.dim; i++)
        {
            p_init[i]=p_new[i];
            q_init[i]=q_new[i];
        }
        
        /* Check no-U-turn criterion. */
        angle1=0.0;
        angle2=0.0;
        for (int i=0; i<in.dim; i++)
        {
            angle1+=p_new[i]*(q_new[i]-q[i]);
            angle2+=p[i]*(q[i]-q_new[i]);
        }
        
        if (angle1<0.0 && angle2<0.0) break;
    }
    
    if (in.verbose) printf("integration steps: %d\n",it);
    
    /* Transcribe values. ---------------------------------------------------------*/
    
    for (int i=0; i<in.dim; i++)
    {
        q[i]=q_new[i];
        p[i]=p_new[i];
    }
    
    /* Clean up. ------------------------------------------------------------------*/
    
    if (output_trajectory) fclose(pfile);
    
    delete[] q_init;
    delete[] p_half;
    delete[] p_init;
    delete[] out;
}

/* Proposal based on the solution of Hamilton's equation. -------------------------*/
void HMC::propose()
{
    /* Draw random prior momenta. */
    for (int i=0; i<in.dim; i++) p[i]=randn(0.0,sqrt(m[i]));
    
    /* Integrate Hamilton's equations. */
    //leap_frog(in.verbose);
}

/* Misfit for Hamiltonian Monte Carlo. --------------------------------------------*/
double HMC::energy()
{
    /* Local variables. */
    double H;
    
    /* Compute energy - model part. */
    H=0.5*in.C;
    
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



/********************************************************************************//**
 * Base class for Monte Carlo sampling.
************************************************************************************/

/* Constructor. -------------------------------------------------------------------*/
mc::mc(Input &in)
{
    verbose = false;
    
    /* General setup. -------------------------------------------------------------*/
    nt = in.hmc_nt;
    dt = in.hmc_dt;
    Nq = in.dim;
    iterations = in.n_samples;
    gamma = in.hmc_gamma;
    m_reg = in.hmc_reg;
    
    /* Set model parameters (mean, standard deviation, scale). --------------------*/
    for (int i=0; i<Nq; i++)
    {
        q.mean_q[i] = in.mean_q[i];
        q.sigma_q[i] = in.sigma_q[i];
        q.q[i] = in.mean_q[i];
        q.scale[i] = in.scale_q[i];
        
        q_new.mean_q[i] = in.mean_q[i];
        q_new.sigma_q[i] = in.sigma_q[i];
        q_new.q[i] = in.mean_q[i];
        q_new.scale[i] = in.scale_q[i];
    }
    
    /* Initialise random number generator. ----------------------------------------*/
    srand (time(0));
    
    /* Pre-computed matrices and vectors. -----------------------------------------*/
    B=new double[Nq];
    A=new double*[Nq];
    for (int i=0; i<Nq; i++) A[i]=new double[Nq];
    
    C = in.C;
    for (int i=0; i<Nq; i++)
    {
        B[i] = in.B[i];
        for (int j=0; j<Nq; j++) A[i][j] = in.A[i][j];
    }
    
    /* Regularisation / prior. ----------------------------------------------------*/
    for (int i=0; i<Nq; i++) A[i][i]+=A[i][i]+0.5/(q.sigma_q[i]*q.sigma_q[i]);
    
    if (verbose)
    {
        printf("A=\n");
        for (int i=0; i<Nq; i++)
        {
            for (int j=0; j<Nq; j++)
            {
                if (A[i][j]>=0) printf("+%.3e ",A[i][j]);
                else printf("%.3e ",A[i][j]);
            }
            printf("\n");
        }
    }
    
    /* Initialise mass matrix. ----------------------------------------------------*/
    m=new double[Nq];
    for (int i=0; i<Nq; i++) m[i]=gamma*(A[i][i]+m_reg);
    
    /* Initialise models and momenta. ---------------------------------------------*/
    p=new double[Nq];
    p_new=new double[Nq];
    
    for (int i=0; i<Nq; i++)
    {
        p_new[i]=randn(0.0,sqrt(m[i]));
        q_new.q[i]=randn(q.mean_q[i],q.sigma_q[i]);
        q.q[i]=q_new.q[i];
    }
}

/* Destructor. --------------------------------------------------------------------*/
mc::~mc()
{
    if (p) delete[] p;
    if (p_new) delete[] p_new;
    if (m) delete[] m;
    if (B) delete[] B;
    if (A)
    {
        for (int i=0; i<Nq; i++) if (A[i]) delete[] A[i];
        delete[] A;
    }
}

/* Debug. -------------------------------------------------------------------------*/
void mc::print()
{
    printf("%d %d %d %lg %lg %lg\n", Nq, iterations, nt, dt, gamma, m_reg);
    
    for (int i=0; i<Nq; i++)
    {
        printf("%lg %lg %lg\n",q.mean_q[i], q.sigma_q[i], q.q[i]);
    }
    
    
    
    printf("\n");
}


/********************************************************************************//**
 * Proposals.
************************************************************************************/

/* Propose test model based on prior. ---------------------------------------------*/
void mc::propose_metropolis()
{
    for (int i=0; i<Nq; i++) q_new.q[i]=randn(q.mean_q[i],q.sigma_q[i]);
}

/* Proposal based on the solution of Hamilton's equation. -------------------------*/
void mc::propose_hamilton()
{
    /* Draw random prior momenta. */
    for (int i=0; i<Nq; i++) p[i]=randn(0.0,sqrt(m[i]));
    
    /* Integrate Hamilton's equations. */
    leap_frog(verbose);
}


/********************************************************************************//**
 * Misfits.
************************************************************************************/

/* Misfit for Metropolis Hastings. ------------------------------------------------*/
double mc::chi()
{
    /* Local variables. */
    double likelihood;
    double a;
    
    /* Compute likelihood. */
    likelihood=0.5*C;
    for (int i=0; i<Nq; i++)
    {
        likelihood+=(q_new.q[i]-q.mean_q[i])*B[i];
        for (int j=0; j<Nq; j++)
        {
            a=A[i][j];
            if (j==i) a-=1.0/(q.sigma_q[i]*q.sigma_q[i]); /* Subtract prior to produce likelihood. */
            likelihood+=0.5*(q_new.q[i]-q.mean_q[i])*(q_new.q[j]-q.mean_q[j])*a;
        }
    }
    
    return likelihood;
}

/* Misfit for Hamiltonian Monte Carlo. --------------------------------------------*/
double mc::energy()
{
    /* Local variables. */
    double H;
    
    /* Compute energy - model part. */

    H=0.5*C;
    
    for (int i=0; i<Nq; i++)
    {
        H+=(q_new.q[i]-q.mean_q[i])*B[i];
        for (int j=0; j<Nq; j++)
        {
            H+=0.5*(q_new.q[i]-q.mean_q[i])*(q_new.q[j]-q.mean_q[j])*A[i][j];
        }
    }
    
    /* Compute energy - momentum part. */
    for (int i=0; i<Nq; i++) H+=0.5*p_new[i]*p_new[i]/m[i];
    
    return H;
}

/* Potential energy, i.e. posterior. ---------------------------------------------*/
double mc::potential_energy()
{
    /* Local variables. */
    double U;
    
    /* Compute energy - model part. */
    
    U=0.5*C;
    
    for (int i=0; i<Nq; i++)
    {
        U+=(q_new.q[i]-q.mean_q[i])*B[i];
        for (int j=0; j<Nq; j++)
        {
            U+=0.5*(q_new.q[i]-q.mean_q[i])*(q_new.q[j]-q.mean_q[j])*A[i][j];
        }
    }
    
    return U;
}


/*=================================================================================*/
/* Right-hand sides of Hamilton equations. ----------------------------------------*/
/*=================================================================================*/

void mc::dHdp(double *p_in, double *rhs)
{
    for (int i=0; i<Nq; i++) rhs[i]=p_in[i]/m[i];
}

void mc::dHdq(parameters &q_in, double *rhs)
{
    /* March through components. */
    for (int k=0; k<Nq; k++)
    {
        rhs[k]=B[k];
        for (int i=0; i<Nq; i++) rhs[k]+=(q_in.q[i]-q_in.mean_q[i])*A[i][k];
    }
}

/********************************************************************************//**
 * Miscellaneous.
************************************************************************************/

/* Write a sample to an open file. ------------------------------------------------*/
void mc::write_sample(FILE *pfile, double misfit, int iteration, int multiplicity)
{
    if (iteration==0) fprintf(pfile,"%d %d\n",Nq,iterations+1);
    
    for (int i=0; i<Nq; i++) fprintf(pfile,"%lg ",q.scale[i]*q.q[i]);
    fprintf(pfile,"%lg %d\n",misfit,multiplicity);
}

/* Leap-frog integration of Hamilton's equations. ---------------------------------*/
void mc::leap_frog(bool verbose)
{
    /* Local variables and setup. -------------------------------------------------*/
    
    int it;
    
    double *p_half, *p_init, *out;
    double angle1, angle2;
    parameters q_init;
    
    out=new double[Nq];
    p_half=new double[Nq];
    p_init=new double[Nq];
    
    FILE *pfile;
    if (verbose) pfile=fopen("OUTPUT/trajectory.txt","a");
    
    /* Set initial values. --------------------------------------------------------*/
    
    q_init=q;
    q_new=q;
    for (int i=0; i<Nq; i++) p_init[i]=p[i];
    
    /* March forward. -------------------------------------------------------------*/
    
    if (verbose)
    {
        //fprintf(pfile,"%d %d\n",2*Nq,nt);
        for (int i=0; i<2*Nq; i++) fprintf(pfile,"0.0 ");
        fprintf(pfile,"\n");
    }
    
    dHdq(q_init,out);
    
    for (it=0; it<nt; it++)
    {
        /* Some output. */
        if (verbose)
        {
            for (int i=0; i<Nq; i++) fprintf(pfile,"%lg ",q_init.q[i]);
            for (int i=0; i<Nq; i++) fprintf(pfile,"%lg ",p_init[i]);
            fprintf(pfile,"\n");
        }
        
        /* First half step in momentum. */
        for (int i=0; i<Nq; i++)
        {
            p_half[i]=p_init[i]-0.5*dt*out[i];
        }
        
        /* Full step in position. */
        dHdp(p_half,out);
        for (int i=0; i<Nq; i++)
        {
            q_new.q[i]=q_init.q[i]+dt*out[i];
        }
        
        /* Second half step in momentum. */
        dHdq(q_new,out);
        for (int i=0; i<Nq; i++)
        {
            p_new[i]=p_half[i]-0.5*dt*out[i];
        }
        
        /* Update position and momentum. */
        for (int i=0; i<Nq; i++)
        {
            p_init[i]=p_new[i];
            q_init.q[i]=q_new.q[i];
        }
        
        /* Check no-U-turn criterion. */
        angle1=0.0;
        angle2=0.0;
        for (int i=0; i<Nq; i++)
        {
            angle1+=p_new[i]*(q_new.q[i]-q.q[i]);
            angle2+=p[i]*(q.q[i]-q_new.q[i]);
        }
        
        if (angle1<0.0 && angle2<0.0) break;
    }
    
    if (verbose) printf("integration steps: %d\n",it);
    
    /* Clean up. ------------------------------------------------------------------*/
    
    if (verbose) fclose(pfile);
    
    delete[] p_half;
    delete[] p_init;
    delete[] out;
}




