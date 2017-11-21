import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('xtick', labelsize=20) 
mpl.rc('ytick', labelsize=20) 

noise=[10.0, 20.0, 50.0, 100.0, 200.0, 1000.0, 5000.0, 10000.0, 50000.0, 100000.0]
mean=[1.335, 0.977, 0.952, 0.975, 1.018, 0.9723, 0.9915, 1.009, 1.0058, 1.0023]
stdv=[-0.74, 0.138, 0.19, 0.1759, 0.155, 0.1747, 0.1427, 0.1689, 0.1851, 0.18]

plt.semilogx(noise, mean, 'ko-',linewidth=1.5)
plt.semilogx(noise, stdv, 'ko--',linewidth=1.5)
plt.grid()
plt.xlabel('number of samples',fontsize=20)
plt.ylabel('posterior mean / standard deviation',fontsize=20)
plt.savefig('OUTPUT/convergence.png')
plt.close()
#plt.show()
