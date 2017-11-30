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
    clock_t start=clock();
    
    /* Local variables. ---------------------------------------------------------------*/
    HMC hmc(argv[1]);
    
    /* Initial values. ----------------------------------------------------------------*/
    for (int i=0; i<hmc.in.dim; i++) hmc.q[i]=hmc.in.mean_q[i];
    for (int i=0; i<hmc.in.dim; i++) hmc.p[i]=0.0;
    
    hmc.p[0]=1.0;
    hmc.q[1]=0.0;
    
    /* Integrate. ---------------------------------------------------------------------*/
    hmc.leap_frog();
    
    printf("elapsed time: %f\n",(double)(clock()-start)/CLOCKS_PER_SEC);
}

