import numpy as np
import random as random
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

#============================================================
#- Setup.
#============================================================

#- Number of burn-in samples to be ignored.
nbi=100

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (8, 8),
         'axes.labelsize': 20,
         'axes.titlesize':'x-large',
         'xtick.labelsize': 20,
         'ytick.labelsize': 20}
pylab.rcParams.update(params)

#============================================================
#- Read samples and plot trajectory.
#============================================================

fid=open('OUTPUT/samples.txt')
dummy=fid.read().strip().split()
fid.close()

dimension=int(dummy[0])
iterations=int(dummy[1])-nbi

bin_start=0.2
bin_end=1.7
bin_width=0.05

bins=np.arange(bin_start,bin_end,bin_width)
hist=np.zeros([len(bins),10])

t=np.arange(0.0,20.0,2.0)
B,T=np.meshgrid(t,bins)

for i in range(100000):#iterations):

	for dim in range(10,20):
		x=float(dummy[2+dim+(i+nbi)*(dimension+1)])

		for bin in range(len(bins)):
			if (x>=bins[bin] and x<bins[bin]+bin_width):
				hist[bin,dim-10]+=1.0


plt.pcolor(B,T,hist,cmap='binary')
plt.xlabel('time')
plt.ylabel('s')
plt.title('source time function')
plt.savefig('OUTPUT/stf.png')
plt.close()
#plt.show()

