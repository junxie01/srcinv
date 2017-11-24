#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mc.hpp"
#include "aux.hpp"

/* ./sampling "input_filename.toml" */

int main(int argc, char *argv[])
{
    /* Check input. -------------------------------------------------------------------*/
    //input in(argv[1]);
    Input in;               //// May be later deleted.
    in.read(argv[1]);
    clock_t start=clock();
    bool verbose = false;
    
    /*---------------------------------------------------------------------------------*/
    /* New Hamiltonian Monte Carlo. ---------------------------------------------------*/
    /*---------------------------------------------------------------------------------*/
    
    if (!strcmp(in.method,"HMC_new"))
    {
        /* Local variables. -----------------------------------------------------------*/
        HMC hmc(argv[1]);
        int accepted = 0;
        int multiplicity = 1;
        double H, H_new, U, U_new;
        FILE *pfile;
        pfile=fopen("../output/samples.txt","w");
        
        /* Initial values. ------------------------------------------------------------*/
        H=hmc.energy();
        U=hmc.potential_energy();
        hmc.write_sample(pfile,U,0,1);
        
        /* Random walk. ---------------------------------------------------------------*/
        for (int it=0; it<hmc.in.n_samples; it++)
        {
            /* Make a model proposition and compute misfit/energy. */
            hmc.propose();
            H_new=hmc.energy();
            
            /* Check Metropolis rule. */
            if ((H_new<H) || (exp(H-H_new)>randf(0.0,1.0)))
            {
                hmc.leap_frog();
                U=hmc.potential_energy();
                hmc.write_sample(pfile,U,it+1,multiplicity);
                
                multiplicity=1;
                H=H_new;
                accepted++;
            }
            else multiplicity++;
        }
        
        printf("accepted: %d\n",accepted);
        fclose(pfile);
    }
    
    /*---------------------------------------------------------------------------------*/
    /* Hamiltonian Monte Carlo. -------------------------------------------------------*/
    /*---------------------------------------------------------------------------------*/
    
    if (!strcmp(in.method,"HMC"))
    {
        /* Local variables. -----------------------------------------------------------*/
        mc m(in);
        int accepted = 0;
        int multiplicity = 1;
        double x, x_new, U, U_new;
        FILE *pfile;
        pfile=fopen("../output/samples.txt","w");
        
        /* Initial values. ------------------------------------------------------------*/
        x=m.energy();
        U=m.potential_energy();
        m.write_sample(pfile,x,0,1);

        /* Random walk. ---------------------------------------------------------------*/
        for (int it=0; it<m.iterations; it++)
        {
            /* Make a model proposition and compute misfit/energy. */
            m.propose_hamilton();
            x_new=m.energy();
            U_new=m.potential_energy();
            
            /* Check Metropolis rule. */
            if ((x_new<x) || (exp(x-x_new)>randf(0.0,1.0)))
            {
                m.write_sample(pfile,U,it+1,multiplicity);
                multiplicity=1;
                
                x=x_new;
                U=U_new;
                m.q=m.q_new;
                accepted++;
            }
            else
            {
                multiplicity++;
            }
        }
        
        printf("accepted: %d\n",accepted);
        fclose(pfile);
    }
    
    /*---------------------------------------------------------------------------------*/
    /* Metropolis Hastings. -----------------------------------------------------------*/
    /*---------------------------------------------------------------------------------*/

    else if (!strcmp(in.method,"MH"))
    {
        /* Local variables. -----------------------------------------------------------*/
        mc m(in);
        int accepted=0;
        double x, x_new, U, U_new;
        FILE *pfile;
        pfile=fopen("../output/samples.txt","w");

        /* Initial values. ------------------------------------------------------------*/
        x=m.chi();
        U=m.potential_energy();
        m.write_sample(pfile,x,0,1);
        
        /* Random walk. ---------------------------------------------------------------*/
        
        for (int it=0; it<m.iterations; it++)
        {
            /* Make a model proposition and compute misfit/energy. */
            m.propose_metropolis();
            x_new=m.chi();
            U_new=m.potential_energy();
            
            /* Check Metropolis rule. */
            if ((x_new<x) || (exp(x-x_new)>randf(0.0,1.0)))
            {
                x=x_new;
                U=U_new;
                m.q=m.q_new;
                accepted++;
            }
            
            m.write_sample(pfile,U,it+1,1);
        }
        
        printf("accepted: %d\n",accepted);
        fclose(pfile);
    }
    
    printf("elapsed time: %f\n",(double)(clock()-start)/CLOCKS_PER_SEC);
}

