import numpy as np
import random as random
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

def plot_samples(dim_1,dim_2,burn_in_samples=50,N=1000,optimal=False,trajectory=False,statistics=False,convergence=False):
    """
    dim_1, dim_2: phase space dimension
    burn_in_samples: number of burn-in samples to be omitted
    N: first N samples used for statistics
    optimal: determine and write optimal model among the samples
    trajectory: plot model-space trajectory
    statistics: compute statistics (standard deviations and covariances)
    convergence: plot convergence towards the final distribution
    """
    #============================================================
    #- Setup.
    #============================================================

    #- Number of burn-in samples to be ignored.
    nbi=burn_in_samples

    #- Incremental displacement for duplicate points.
    epsilon_1=0.0003
    epsilon_2=0.0003
    
    #- Names of parameters.
    names=["Mxx","Mxy","Mxz","Myy","Myz","Mzz","lon","lat","z","t0","s0","s1","s2","s3","s4","s5","s6","s7","s8","s9"]

    params = {'legend.fontsize': 'x-large',
          'figure.figsize': (8, 8),
         'axes.labelsize': 20,
         'axes.titlesize':'x-large',
         'xtick.labelsize': 20,
         'ytick.labelsize': 20}
    pylab.rcParams.update(params)

    #- Write parameter names.
    print '========= ', names[dim_1], names[dim_2], '============'

    #============================================================
    #- Read samples and plot trajectory.
    #============================================================

    fid=open('OUTPUT/samples.txt','r')
    dummy=fid.read().strip().split()
    fid.close()

    dimension=int(dummy[0])
    iterations=int(dummy[1])-nbi

    x=np.zeros(iterations)
    y=np.zeros(iterations)

    x_plot=np.zeros(iterations)
    y_plot=np.zeros(iterations)

    q_opt=np.zeros(dimension)
    chi=1.0e100

    for i in range(iterations):

        x[i]=float(dummy[2+dim_1+(i+nbi)*(dimension+1)])
        y[i]=float(dummy[2+dim_2+(i+nbi)*(dimension+1)])
        x_plot[i]=x[i];
        y_plot[i]=y[i];

        chi_test=float(dummy[2+(dimension)+(i+nbi)*(dimension+1)])
        if (chi_test<chi):
            chi=chi_test
            #print i, ': chi_min=', chi
            for k in range(dimension):
                q_opt[k]=float(dummy[2+k+(i+nbi)*(dimension+1)])
	
        if (i>0 and x[i]==x[i-1] and y[i]==y[i-1]):
            x_plot[i]+=epsilon_1*random.gauss(0.0,1.0)
            y_plot[i]+=epsilon_2*random.gauss(0.0,1.0)

    if trajectory:

        plt.plot(x_plot,y_plot,'k')
        plt.plot(x_plot,y_plot,'ro')
        plt.axis('equal')
        plt.xlabel(names[dim_1])
        plt.ylabel(names[dim_2])
        plt.title('random walk trajectory')
        plt.savefig('OUTPUT/trajectory.png')
        plt.close()

    #============================================================
    #- Optimal model.
    #============================================================

    if optimal:

        print "-- optimal model --"

        fid=open('OUTPUT/optimal_model.txt','w')

        for k in range(dimension):
            print names[k]+"="+str(q_opt[k])
            fid.write(str(q_opt[k])+'\n')

        fid.close()

    #============================================================
    #- Histograms.
    #============================================================

    xlim=np.max(np.abs(x));
    ylim=np.max(np.abs(y));

    plt.hist(x,bins=25,color='k',normed=True)
    #plt.xlim([-xlim,xlim])
    plt.xlabel(names[dim_1])
    plt.ylabel('posterior marginal')
    plt.savefig('OUTPUT/marginal_'+names[dim_1]+'.png')
    plt.close()

    plt.hist(y,bins=25,color='k',normed=True)
    #plt.xlim([-ylim,ylim])
    plt.xlabel(names[dim_2])
    plt.ylabel('posterior marginal')
    plt.savefig('OUTPUT/marginal_'+names[dim_2]+'.png')
    plt.close()

    plt.hist2d(x,y,bins=15,normed=True,cmap='binary')
    #plt.axis('equal')
    plt.xlabel(names[dim_1])
    plt.ylabel(names[dim_2])
    plt.title('2D posterior marginal')
    plt.colorbar()
    plt.savefig('OUTPUT/marginal_2D_'+names[dim_1]+'_'+names[dim_2]+'.png')
    plt.close()

    #============================================================
    #- Assess convergence.
    #============================================================

    if convergence:

        n=range(nbi,iterations,10)

        hist_final,bin=np.histogram(x,bins=20,density=True)
        diff=np.zeros(len(n))

        k=0
        for i in n:
            hist,bin=np.histogram(x[0:i],bins=20,density=True)
            diff[k]=np.sum(np.abs(hist-hist_final))
            k=k+1

        plt.semilogy(n,diff,'k')
        plt.xlabel('samples')
        plt.ylabel('difference to final')
        plt.savefig('OUTPUT/convergence1.png')
        plt.close()


        hist_final,bin=np.histogram(y,bins=20,density=True)
        diff=np.zeros(len(n))

        k=0
        for i in n:
            hist,bin=np.histogram(y[0:i],bins=20,density=True)
            diff[k]=np.sum(np.abs(hist-hist_final))
            k=k+1

        plt.semilogy(n,diff,'k')
        plt.xlabel('samples')
        plt.ylabel('difference to final')
        plt.savefig('OUTPUT/convergence2.png')
        plt.close()

    #============================================================
    #- Make statistics.
    #============================================================

    if statistics:

        print "-- statistics --"

        mean_x_all=np.mean(x)
        mean_y_all=np.mean(y)

        mean_x_N=np.mean(x[:N])
        mean_y_N=np.mean(y[:N])

        cov_xx_all=0.0
        cov_yy_all=0.0
        cov_xy_all=0.0

        cov_xx_N=0.0
        cov_yy_N=0.0
        cov_xy_N=0.0

        for i in range(iterations):
            cov_xx_all+=(mean_x_all-x[i])*(mean_x_all-x[i])
            cov_yy_all+=(mean_y_all-y[i])*(mean_y_all-y[i])
            cov_xy_all+=(mean_x_all-x[i])*(mean_y_all-y[i])

        for i in range(N):
            cov_xx_N+=(mean_x_N-x[i])*(mean_x_N-x[i])
            cov_yy_N+=(mean_y_N-y[i])*(mean_y_N-y[i])
            cov_xy_N+=(mean_x_N-x[i])*(mean_y_N-y[i])

        cov_xx_all=cov_xx_all/(iterations)
        cov_yy_all=cov_yy_all/(iterations)
        cov_xy_all=cov_xy_all/(iterations)

        cov_xx_N=cov_xx_N/N
        cov_yy_N=cov_yy_N/N
        cov_xy_N=cov_xy_N/N

        print '--- all samples ---'
        print 'mean_x=', mean_x_all, 'mean_y=', mean_y_all
        print 'std_xx=', np.sqrt(cov_xx_all), 'std_yy=', np.sqrt(cov_yy_all), 'cov_xy=', cov_xy_all/(np.sqrt(cov_xx_all)*np.sqrt(cov_yy_all))

        print '--- first ', N, ' samples ---'
        print 'mean_x=', mean_x_N, 'mean_y=', mean_y_N
        print 'std_xx=', np.sqrt(cov_xx_N), 'std_yy=', np.sqrt(cov_yy_N), 'cov_xy=', cov_xy_N/(np.sqrt(cov_xx_N)*np.sqrt(cov_yy_N))

