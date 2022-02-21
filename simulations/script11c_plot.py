#### Study the effect of the number of directions

import numpy as np
import pickle as pk
from matplotlib import pyplot as plt
from dipy.io import read_bvals_bvecs
from dipy.data import gradient_table
import dipy.reconst.dki as dki
from AdcAkcFit import voxelfit as vfit


## Input information: number of gradient directions to study
dirs=[3,9,21,30,61]
ndirs = len(dirs)
sigfiles=[]
for nn in range(0,ndirs):
	sigfiles.append('results_script11/dir{}/syn_d20.0_D75.0_bmin100.0_bmax2000.0_Nbval7.sig.icell.SNRinf.bin'.format(dirs[nn]))
msize = [5,8,11,14,17]    # Marker size for the different number of directions
snr = 20                  # SNR values for simulations

##### SIGNALS

## Plot: show signals after averaging over all gradient directions
plt.subplot(2,3,1)
for nn in range(0,ndirs):

	# Load MRI signals and corresponding diffusion MRI protocols
	h = open(sigfiles[nn],'rb')
	data = pk.load(h)
	h.close()
	data = data[0]
	sig = data[0]
	prot = data[1]
	bvals = prot[0]*1000.0   # b-values in s/mm2
	gx = prot[4]             # Gradient directions along x
	gy = prot[5]             # Gradient directions along y
	gz = prot[6]             # Gradient directions along z
	
	# Calculate powder average
	bunique = np.unique(bvals)
	sunique = np.zeros(bunique.size)
	for bb in range(0,bunique.size):
		sunique[bb] = np.mean(sig[bvals==bunique[bb]])
	bunique = np.concatenate(([0.0],bunique))
	sunique = np.concatenate(([1.0],sunique))
	
	plt.plot(bunique,sunique,'o',markersize=msize[nn],markerfacecolor='none',label='{} directions'.format(dirs[nn]))

plt.yticks(ticks=[0.4,0.5,0.6,0.7,0.8,0.9,1.0],fontsize=14)
plt.xticks(ticks=bunique,labels=[0,0.1,0.4,0.7,1.0,1.4,1.7,2.0],fontsize=10)
plt.ylim([0.795,1.01])
plt.grid()
plt.xlabel('$b$   [ms/$\mu$m$^2$]',fontsize=12)
plt.ylabel('Directionally-averaged signal   $s(b)$',fontsize=13)
plt.title('(A) Directional averaging, SNR→∞',fontsize=13,fontweight='bold')
plt.legend(fontsize=12)		


##### DIFFUSIVITY

# Boxplot: show apparent diffusion coefficient calculated with custom-written code on directionally-averaged measurements
hfile = open('./script11b_analyse.adc.bin'.format(snr),'rb')
adcdata = pk.load(hfile) 
hfile.close()
adcdata = np.transpose(adcdata)    # D is defined in um2/ms
plt.subplot(2,3,2)
plt.boxplot(adcdata)
plt.ylabel('$D$  [$\mu$m$^2$/ms]',fontsize=13)
plt.ylim([0.0,0.3])
#plt.yticks(ticks=[-4,-2,0,2,4,6,8,10])
plt.xticks(ticks=[1,2,3,4,5],labels=['3 dir.','9 dir.','21 dir.','30 dir.','61 dir.'],rotation=45,fontsize=12)
plt.title('(B) $D$, Directional averaging',fontsize=13,fontweight='bold')
plt.grid()

# Boxplot: show mean diffusivity from full kurtosis tensor fit performed with DiPy
hfile = open('./script11b_analyse.md.bin'.format(snr),'rb')
mddata = pk.load(hfile) 
hfile.close()
mddata = 1000.0*np.transpose(mddata)    # Convert units of MD to um2/ms, as for D
plt.subplot(2,3,3)
plt.boxplot(mddata)
plt.ylabel('$MD$  [$\mu$m$^2$/ms]',fontsize=13)
plt.ylim([0.0,0.3])
#plt.yticks(ticks=[-4,-2,0,2,4,6,8,10])
plt.xticks(ticks=[1,2,3,4,5],labels=['3 dir.','9 dir.','21 dir.','30 dir.','61 dir.'],rotation=45,fontsize=12)
plt.title('(C) $MD$, Tensor fit',fontsize=13,fontweight='bold')
plt.grid()




##### KURTOSIS

# Boxplot: show apparent kurtosis coefficient calculated with custom-written code on directionally-averaged measurements
hfile = open('./script11b_analyse.akc.bin'.format(snr),'rb')
akcdata = pk.load(hfile) 
hfile.close()
akcdata = np.transpose(akcdata)
plt.subplot(2,3,5)
plt.boxplot(akcdata)
plt.ylabel('$K$',fontsize=13)
plt.ylim([-5.2,10.5])
plt.yticks(ticks=[-4,-2,0,2,4,6,8,10])
plt.xticks(ticks=[1,2,3,4,5],labels=['3 dir.','9 dir.','21 dir.','30 dir.','61 dir.'],rotation=45,fontsize=12)
plt.title('(D) $K$, Directional averaging',fontsize=13,fontweight='bold')
plt.grid()

# Boxplot: show mean kurtosis from full kurtosis tensor fit performed with DiPy
hfile = open('./script11b_analyse.mk.bin'.format(snr),'rb')
mkdata = pk.load(hfile) 
hfile.close()
mkdata = np.transpose(mkdata)
plt.subplot(2,3,6)
plt.boxplot(mkdata)
plt.ylabel('$MK$',fontsize=13)
plt.ylim([-5.2,10.5])
plt.yticks(ticks=[-4,-2,0,2,4,6,8,10])
plt.xticks(ticks=[1,2,3,4,5],labels=['3 dir.','9 dir.','21 dir.','30 dir.','61 dir.'],rotation=45,fontsize=12)
plt.title('(E) $MK$, Tensor fit',fontsize=13,fontweight='bold')
plt.grid()


# Set up figure
plt.subplots_adjust(left=0.06, bottom=0.09, right=0.98, top=0.94, wspace=0.25, hspace=0.35)
fig = plt.gcf()
fig.set_size_inches([14.0,8.0])

# Show and save figure
plt.savefig('script11c_plot.900dpi.png',dpi=900)
plt.show()


