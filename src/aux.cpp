
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include "aux.hpp"
#include "toml.h"

const double PI = 3.14159265358979323846264338327;

/*=================================================================================*/
/* Input --------------------------------------------------------------------------*/
/*=================================================================================*/

/* Constructor. -------------------------------------------------------------------*/
input::input(const char *input_filename)
{
    /* Local variables. -----------------------------------------------------------*/
    FILE *fid;
    char str[1000];
    toml::Value* x;
    std::vector<double> y;
    toml::ParseResult pr = toml::parseFile(input_filename);
    
    /* General setup and model parameters. ----------------------------------------*/
    
    /* Model space dimension. */
    x = pr.value.find("dim");
    dim = x->as<int>();
    
    /* Method. */
    x = pr.value.find("method");
    std::string s = x->as<std::string>();
    strcpy(method,s.c_str());
    
    /* Allocate memory before filling values. */
    mean_q = new double[dim];
    sigma_q = new double[dim];
    scale_q = new double[dim];
    
    /* Moment tensor components, prior means. */
    x = pr.value.find("M");
    y = x->as<std::vector<double>>();
    for (int i=0; i<6; i++) mean_q[i] = y[i];

    /* Moment tensor components, prior standard deviations. */
    x = pr.value.find("sigma_M");
    y = x->as<std::vector<double>>();
    for (int i=0; i<6; i++) sigma_q[i] = y[i];

    /* Moment tensor components scalings. */
    x = pr.value.find("scale_M");
    y = x->as<std::vector<double>>();
    for (int i=0; i<6; i++) scale_q[i] = y[i];

    if (dim > 6)
    {
        /* Source location. */
        mean_q[6] = pr.value.find("lon")->as<double>();
        mean_q[7] = pr.value.find("lat")->as<double>();
        mean_q[8] = pr.value.find("depth")->as<double>();
        
        sigma_q[6] = pr.value.find("sigma_lon")->as<double>();
        sigma_q[7] = pr.value.find("sigma_lat")->as<double>();
        sigma_q[8] = pr.value.find("sigma_depth")->as<double>();
        
        scale_q[6] = pr.value.find("scale_lon")->as<double>();
        scale_q[7] = pr.value.find("scale_lat")->as<double>();
        scale_q[8] = pr.value.find("scale_depth")->as<double>();

    }
    
    /* Origin time. */
    if (dim == 10)
    {
        mean_q[9] = pr.value.find("t")->as<double>();
        sigma_q[9] = pr.value.find("sigma_t")->as<double>();
        scale_q[9] = pr.value.find("scale_t")->as<double>();
    }
    
    /* Monte Carlo-specific parameters. -------------------------------------------*/
    n_samples = pr.value.find("n_samples")->as<int>();
    
    /* HMC tuning parameters. -----------------------------------------------------*/
    hmc_dt = pr.value.find("hmc_dt")->as<double>();
    hmc_nt = pr.value.find("hmc_nt")->as<int>();
    hmc_gamma = pr.value.find("hmc_gamma")->as<double>();
    hmc_reg = pr.value.find("hmc_reg")->as<double>();
    
    /* Matrix, vector and constant. -----------------------------------------------*/
    
    /* Initialise. */
    B=new double[dim];
    A=new double*[dim];
    for (int i=0; i<dim; i++) A[i]=new double[dim];
    
    /* Read C. */
    fid=fopen("../input/c.txt","r");
    fscanf(fid,"%lg",&C);
    fclose(fid);
    
    /* Read B[i]. */
    fid=fopen("../input/b.txt","r");
    for (int i=0; i<dim; i++)
    {
        fscanf(fid,"%lg",&(B[i]));
        B[i]=B[i]*scale_q[i];
    }
    fclose(fid);
    
    /* Read A[i][j]. */
    fid=fopen("../input/A.txt","r");
    for (int i=0; i<dim; i++)
    {
        for (int j=0; j<dim; j++)
        {
            fscanf(fid,"%lg",&(A[i][j]));
            A[i][j]=A[i][j]*scale_q[i]*scale_q[j];
        }
        fgets(str,1000,fid);
    }
    fclose(fid);

}

/* Destructor. --------------------------------------------------------------------*/
input::~input()
{
    if (mean_q) delete[] mean_q;
    if (sigma_q) delete[] sigma_q;
    if (scale_q) delete[] scale_q;
    
    if (B) delete[] B;
    if (A)
    {
        for (int i=0; i<dim; i++) if (A[i]) delete[] A[i];
        delete[] A;
    }
}


/*=================================================================================*/
/* Modelling parameters class. ----------------------------------------------------*/
/*=================================================================================*/

/* Constructor. -------------------------------------------------------------------*/
parameters::parameters()
{

}

/* Copy constructor. --------------------------------------------------------------*/

parameters::parameters(const parameters &a)
{
    for (int i=0; i<20; i++)
    {
        q[i]=a.q[i];
        sigma_q[i]=a.sigma_q[i];
        mean_q[i]=a.mean_q[i];
    }
}

/* Assignment operator. -----------------------------------------------------------*/

parameters & parameters::operator=(const parameters &a)
{
    /* Check for self-assignment. */
    if (this == &a) return *this;
    
    for (int i=0; i<20; i++)
    {
        q[i]=a.q[i];
        sigma_q[i]=a.sigma_q[i];
        mean_q[i]=a.mean_q[i];
    }
    
    return *this;
}

/* Destructor. --------------------------------------------------------------------*/
parameters::~parameters()
{
    
}

/* Read input from file. ----------------------------------------------------------*/
void parameters::read_input(void)
{
    FILE *fid;
    char str[1000];
    
    fid=fopen("../input/parameters_source.txt","r");
    
    /* Source coordinates (x, y, z). */
    fgets(str,1000,fid);
    for (int i=6; i<9; i++) fscanf(fid,"%lg",&(mean_q[i]));
    fgets(str,1000,fid);
    for (int i=6; i<9; i++) fscanf(fid,"%lg",&(sigma_q[i]));
    fgets(str,1000,fid);

    /* Moment tensor components. */
    fgets(str,1000,fid);
    for (int i=0; i<6; i++) fscanf(fid,"%lg",&(mean_q[i]));
    fgets(str,1000,fid);
    for (int i=0; i<6; i++) fscanf(fid,"%lg",&(sigma_q[i]));
    fgets(str,1000,fid);
    
    /* Origin time. */
    fgets(str,1000,fid);
    fscanf(fid,"%lg",&(mean_q[9]));
    fgets(str,1000,fid);
    fscanf(fid,"%lg",&(sigma_q[9]));
    fgets(str,1000,fid);

    /* Source time function coefficients. */
    fgets(str,1000,fid);
    for (int i=10; i<20; i++) fscanf(fid,"%lg",&(mean_q[i]));
    fgets(str,1000,fid);
    for (int i=10; i<20; i++) fscanf(fid,"%lg",&(sigma_q[i]));
    fgets(str,1000,fid);
    
    /* Assign values equal to prior means. */
    for (int i=0; i<20; i++) q[i]=mean_q[i];
    
    fclose(fid);
    
    /* Read scaling parameters. */
    fid=fopen("../input/scale.txt","r");
    for (int i=0; i<20; i++) fscanf(fid,"%lg",&(scale[i]));
    fclose(fid);
    
}


/* Addition. ----------------------------------------------------------------------*/

parameters operator+(const parameters &a, const parameters &b)
{
    parameters c;
    
    for (int i=0; i<20; i++)
    {
        c.q[i]=a.q[i]+b.q[i];
        c.mean_q[i]=a.mean_q[i];
        c.sigma_q[i]=a.sigma_q[i];
    }
    
    return c;
}

/* Subtraction. -------------------------------------------------------------------*/

parameters operator-(const parameters &a, const parameters &b)
{
    parameters c;
    
    for (int i=0; i<20; i++)
    {
        c.q[i]=a.q[i]-b.q[i];
        c.mean_q[i]=a.mean_q[i];
        c.sigma_q[i]=a.sigma_q[i];
    }
    
    return c;
}


/*=================================================================================*/
/* Little helpers. ----------------------------------------------------------------*/
/*=================================================================================*/

/* Uniformly distributed, double-valued random numbers. ---------------------------*/
double randf(double min, double max)
{
    return (max-min)*(double)rand()/RAND_MAX+min;
}

/* Normally distributed, double-valued random numbers. ----------------------------*/
/* This function implements the Box-Muller transform to obtain a pair of
 normally distributed random numbers with a given mean and standard deviation. */
void randn(double mean, double stdv, double *x1, double *x2)
{
    double z1=(double)rand()/RAND_MAX;
    double z2=(double)rand()/RAND_MAX;
    
    *x1=sqrt(-2.0*log(z1))*cos(2.0*PI*z2);
    *x2=sqrt(-2.0*log(z1))*sin(2.0*PI*z2);
    
    *x1=stdv*(*x1)+mean;
    *x2=stdv*(*x2)+mean;
}

double randn(double mean, double stdv)
{
    double x;
    
    double z1=(double)rand()/RAND_MAX;
    double z2=(double)rand()/RAND_MAX;
    
    x=sqrt(-2.0*log(z1))*cos(2.0*PI*z2);
    x=stdv*x+mean;
    
    return x;
}


