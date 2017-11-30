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
    Input in;
    in.read(argv[1]);
    clock_t start=clock();
    
    /*---------------------------------------------------------------------------------*/
    /* Hamiltonian Monte Carlo (Fichtner). --------------------------------------------*/
    /*---------------------------------------------------------------------------------*/
    
    if (!strcmp(in.method,"HMC_F"))
    {
        /* Local variables. -----------------------------------------------------------*/
        HMC hmc(argv[1]);
        int accepted = 0;
        int multiplicity = 1;
        double H, H_new;
        FILE *pfile;
        pfile=fopen("../output/samples.txt","w");
        
        /* Initial values. ------------------------------------------------------------*/
        H=hmc.energy();
        hmc.write_sample(pfile,hmc.potential_energy(),0,1);
        
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
                hmc.write_sample(pfile,hmc.potential_energy(),it+1,multiplicity);
                
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
    /* Hamiltonian Monte Carlo (Neal). ------------------------------------------------*/
    /*---------------------------------------------------------------------------------*/
    
    if (!strcmp(in.method,"HMC_N"))
    {
        /* Local variables. -----------------------------------------------------------*/
        HMC hmc(argv[1]);
        int accepted = 0;
        int multiplicity = 1;
        double H, H_new;
        FILE *pfile;
        pfile=fopen("../output/samples.txt","w");
        
        double* q_save=new double[hmc.in.dim];
        
        /* Initial values. ------------------------------------------------------------*/
        H=hmc.energy();
        hmc.write_sample(pfile,hmc.potential_energy(),0,1);
        
        /* Random walk. ---------------------------------------------------------------*/
        for (int it=0; it<hmc.in.n_samples; it++)
        {
        
            /* Make a model proposition and compute misfit/energy. */
            hmc.propose();
            H=hmc.energy();
            for (int i=0; i<hmc.in.dim; i++) q_save[i]=hmc.q[i];
            hmc.leap_frog();
            H_new=hmc.energy();
        
            if ((H_new<H) || (exp(H-H_new)>randf(0.0,1.0)))
            {
                    hmc.write_sample(pfile,hmc.potential_energy(),it+1,multiplicity);
                    multiplicity=1;
                    accepted++;
            }
            else
            {
                for (int i=0; i<hmc.in.dim; i++) hmc.q[i]=q_save[i];
                multiplicity++;
            }
        }
        
        delete[] q_save;
        printf("accepted: %d\n",accepted);
        fclose(pfile);
    }


    /*---------------------------------------------------------------------------------*/
    /* Metropolis Hastings. -----------------------------------------------------------*/
    /*---------------------------------------------------------------------------------*/
    
    else if (!strcmp(in.method,"MH"))
    {
        /* Local variables. -----------------------------------------------------------*/
        MH mh(argv[1]);
        int accepted=0;
        int multiplicity=1;
        double L, L_new, U;
        FILE *pfile;
        pfile=fopen("../output/samples.txt","w");
        
        /* Initial values. ------------------------------------------------------------*/
        L=mh.likelihood();
        U=mh.misfit();
        mh.write_sample(pfile,U,0,1);
        
        /* Random walk. ---------------------------------------------------------------*/
        
        for (int it=0; it<mh.in.n_samples; it++)
        {
            /* Make a model proposition and compute misfit/energy. */
            mh.propose();
            L_new=mh.likelihood();
            
            /* Check Metropolis rule. */
            if ((L_new<L) || (exp(L-L_new)>randf(0.0,1.0)))
            {
                for (int i=0; i<mh.in.dim; i++) mh.q[i]=mh.q_new[i];
                U=mh.misfit();
                mh.write_sample(pfile,U,it+1,multiplicity);
                multiplicity=1;
                L=L_new;
                accepted++;
            }
            else multiplicity++;
        }
        
        printf("accepted: %d\n",accepted);
        fclose(pfile);
    }

    printf("elapsed time: %f\n",(double)(clock()-start)/CLOCKS_PER_SEC);
}

