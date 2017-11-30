import numpy as np
import random as random
import matplotlib.pyplot as plt

def plot_trajectory(dim_1,dim_2,plot=True):

    #============================================================
    #- Setup.
    #============================================================

    #- Determine dimension of phase space.

    fid=open('../output/trajectory.txt')
    s=fid.readline()
    dimension=len(s.split(' '))-1
    fid.close()
    
    #- Names of parameters.
    
    names=dimension*[None]
    names_parameters=["Mxx","Mxy","Mxz","Myy","Myz","Mzz","lon","lat","z","t0"]
    names_momenta=["p_Mxx","p_Mxy","p_Mxz","p_Myy","p_Myz","p_Mzz","p_lon","p_lat","p_z","p_t0"]
    
    for i in range(dimension/2):
        names[i]=names_parameters[i]
        names[i+dimension/2]=names_momenta[i]
    
    #============================================================
    #- Read samples and plot trajectory..
    #============================================================

    #- Plot trajectory. A new sub-trajectory (equal to one HMC iteration) begins with a line of zeros.

    fid=open('../output/trajectory.txt')

    dummy=fid.read().strip().split()
    iterations=len(dummy)/dimension
    fid.close()

    x=np.zeros(iterations)
    y=np.zeros(iterations)

    for i in range(iterations):

        x[i]=float(dummy[dim_1+i*dimension])
        y[i]=float(dummy[dim_2+i*dimension])

    plt.plot(x,y,c=np.random.rand(3,1))
    plt.plot(x,y,'ko')

    plt.xlabel(names[dim_1])
    plt.ylabel(names[dim_2])

    if (plot):
        plt.show()
    else:
        plt.savefig('../output/Hamiltonian_trajectory.png')
        plt.close()
