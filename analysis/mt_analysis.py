from obspy.imaging.beachball import beachball
import numpy as np

#- Input. ---------------------------------------

if True:

	#- Posterior
	Mxx=-7.2924e+15
	Myy=-3.94173e+16
	Mzz=1.54331e+16
	Mxy=-4.02422e+14
	Mxz=-3.92897e+16
	Myz=1.54331e+16

if False:

	#- Prior
	Mxx=-0.778e16
	Myy=-4.279e16
	Mzz=5.0577e16
	Mxy=-0.176e16
	Mxz=1.465e16
	Myz=4.344e16

M=np.array([[Mxx, Mxy, Mxz],[Mxy, Myy, Myz],[Mxz, Myz, Mzz]])

#- Analysis. ------------------------------------

M_iso=np.eye(3)*np.trace(M)/3.0
M_dev=M-M_iso

M0=np.sqrt(Mxx**2+Myy**2+Mzz**2+2*Mxy**2+2*Mxz**2+2*Myz**2)
Mw=(np.log10(M0)-9.05)/1.5

eig=np.linalg.eig(M_dev)[0]

imax=np.argwhere(np.abs(eig)==np.max(np.abs(eig)))
imin=np.argwhere(np.abs(eig)==np.min(np.abs(eig)))

F=-float(eig[imin]/eig[imax])

print "seismic moment, M0="+str(M0)+" Nm"
print "moment magnitude, Mw="+str(Mw)
print "dc percentage: "+str((1.0-2.0*np.abs(F))*100.0)+" %"
print "relative isotropic component: "+str(100.0*np.trace(M)/(M0*np.sqrt(3.0)))+" %"

#- Plot beachball. ------------------------------

mt = [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
beachball(mt, size=200, linewidth=2, facecolor='b')
