import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('xtick', labelsize=20) 
mpl.rc('ytick', labelsize=20) 

noise=[1.0, 0.75, 0.5, 0.25, 0.1, 0.05]
metropolis_noise=[0.68199, 0.53335, 0.31063, 0.06988, 0.00559, 0.00155]
hamilton_noise=[0.40286, 0.4002, 0.39869, 0.40451, 0.40036, 0.4025]

stdv=[0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]
metropolis_stdv=[0.3192, 0.14772, 0.0758, 0.04389, 0.02606, 0.01886, 0.01237, 0.00862, 0.00635]
hamilton_stdv=[0.402, 0.406, 0.401, 0.400, 0.401, 0.402, 0.400, 0.402, 0.400]

plt.semilogy(noise, metropolis_noise, 'ko-',linewidth=1.5)
plt.semilogy(noise, hamilton_noise, 'ko--',linewidth=1.5)
plt.grid()
plt.xlabel('sigma relative to max amplitude',fontsize=20)
plt.ylabel('acceptance rate',fontsize=20)
plt.savefig('OUTPUT/acceptance_noise.png')
plt.close()
#plt.show()

plt.semilogy(stdv, metropolis_stdv, 'ko-',linewidth=1.5)
plt.semilogy(stdv, hamilton_stdv, 'ko--',linewidth=1.5)
plt.grid()
plt.xlabel('stdv',fontsize=20)
plt.ylabel('acceptance rate',fontsize=20)
plt.savefig('OUTPUT/acceptance_stdv.png')
plt.close()
#plt.show()