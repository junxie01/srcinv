import numpy as np
import random as random
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def plot_source(burn_in_samples=50, N=1000):

    #============================================================
    #- Setup.
    #============================================================

    #- Number of burn-in samples to be ignored.
    nbi=burn_in_samples

    #- Names of parameters.
    names=["Mxx","Mxy","Mxz","Myy","Myz","Mzz","lon","lat","z","t0","iso","M0"]

    params = {'legend.fontsize': 'x-large',
          'figure.figsize': (8, 8),
         'axes.labelsize': 20,
         'axes.titlesize':'x-large',
         'xtick.labelsize': 20,
         'ytick.labelsize': 20}
    pylab.rcParams.update(params)

    #================================================================
    #- Read samples and plot trajectory.
    #================================================================

    fid=open('OUTPUT/samples.txt','r')
    dummy=fid.read().strip().split()
    fid.close()

    dimension=int(dummy[0])
    iterations=int(dummy[1])-nbi

    for dim in np.arange(12):
    
        #- Fill array with values. ----------------------------------

        #- Regular parameters.
        if dim<10:
            x=dummy[(2+dim+(nbi)*(dimension+1)):-1:dimension+1]
            x=np.array(x,dtype=float)
        #- Isotropic component.
        if dim==10:
            xx=dummy[(2+0+(nbi)*(dimension+1)):-1:dimension+1]
            yy=dummy[(2+3+(nbi)*(dimension+1)):-1:dimension+1]
            zz=dummy[(2+5+(nbi)*(dimension+1)):-1:dimension+1]
            x=np.array(xx,dtype=float)+np.array(yy,dtype=float)+np.array(zz,dtype=float)
        if dim==11:
            xx=np.array(dummy[(2+0+(nbi)*(dimension+1)):-1:dimension+1],dtype=float)
            xy=np.array(dummy[(2+1+(nbi)*(dimension+1)):-1:dimension+1],dtype=float)
            xz=np.array(dummy[(2+2+(nbi)*(dimension+1)):-1:dimension+1],dtype=float)
            yy=np.array(dummy[(2+3+(nbi)*(dimension+1)):-1:dimension+1],dtype=float)
            yz=np.array(dummy[(2+4+(nbi)*(dimension+1)):-1:dimension+1],dtype=float)
            zz=np.array(dummy[(2+5+(nbi)*(dimension+1)):-1:dimension+1],dtype=float)
            x=np.sqrt(xx**2+yy**2+zz**2+2*xy**2+2*xz**2+2*yz**2)


        #- Normal distribution statistics. --------------------------

        #- All samples.
        mean_all=np.mean(x)
        std_all=np.std(x)

        #- First N samples.
        mean_N=np.mean(x[0:N])
        std_N=np.std(x[0:N])

        print names[dim]+': mean='+str(mean_all)+', std='+str(std_all)+' ('+str(mean_N)+', '+str(std_N)+')'

        #- Plotting. ------------------------------------------------

        n,bins,patches=plt.hist(x,bins=30,color='k',normed=True)
        y_all=mlab.normpdf(bins,mean_all,std_all)
        y_N=mlab.normpdf(bins,mean_N,std_N)
        plt.plot(bins,y_all,'b',linewidth=4.0)
        plt.plot(bins,y_N,'r',linewidth=4.0)
        plt.xlabel(names[dim])

        if dim<6:
            plt.xlim([mean_all-0.8e17, mean_all+0.8e17])
        if (dim==6):
            plt.xlim([mean_all-0.054/0.85, mean_all+0.054/0.85])
        if (dim==7):
            plt.xlim([mean_all-0.054, mean_all+0.054])
        if (dim==8):
            plt.xlim([mean_all-6000.0, mean_all+6000.0])


        plt.ylabel('posterior marginal')
        plt.savefig('OUTPUT/'+names[dim]+'.png')
        plt.close()
