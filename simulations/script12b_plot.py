#### Study the effect of the number of b-values

import numpy as np
import pickle as pk
from matplotlib import pyplot as plt
from dipy.io import read_bvals_bvecs
from dipy.data import gradient_table
import dipy.reconst.dki as dki
from AdcAkcFit import voxelfit as vfit


## Input information: number of gradient directions to study
bmax_array = [1000,2000]
plt_count_title = ['(A)','(B)','(C)','(D)']
snrmax = 100    # Max SNR
snrmin = 90     # Min SNR
subpar = [1,2,3,6,9]   # Parameter to be used to downsample b-values and signals
nsub = len(subpar)   # Number of downsampling experiments
nnoise = 1000    # Number of independent noise instantations
myseed = 100287  # Seed for reproducibility
LW = 3   # Line width


### Loop through maximum number of b-values
plt_count = 1
for bb in range(0,len(bmax_array)):


	## Get maximum b-value
	bmax = bmax_array[bb]

	## ADC and AKC comparison: different protocols with different number of b-values
	sigfile = 'results_script12/bmax{}/syn_d20.0_D75.0_bmin100.0_bmax{}.0_Nbval19.sig.icell.SNRinf.bin'.format(bmax,bmax)
	bvalfile = 'results_script12/bmax{}/syn_d20.0_D75.0_bmin100.0_bmax{}.0_Nbval19.bval'.format(bmax,bmax)
	np.random.seed(myseed)
	adcdata = np.zeros((nsub,nnoise))    # Array storing custom ADC for different noise instantiations and number of b-values
	akcdata = np.zeros((nsub,nnoise))    # Array storing custom AKC for different noise instantiations and number of b-values
	myxlabels = []
	for nn in range(0,nsub):

		# Load MRI signals and corresponding diffusion MRI protocols
		h = open(sigfile,'rb')
		data = pk.load(h)
		h.close()
		data = data[0]
		sig = data[0]
		bvals = np.loadtxt(bvalfile)
		
		
		# Downsample signal and b-values
		sig = sig[0:sig.size:subpar[nn]]
		bvals = bvals[0:bvals.size:subpar[nn]]
		Ndown = sig.size
		myxlabels.append('{} b-values'.format(Ndown))
		
		# Loop through random noise instantiations
		snr_rndvals = np.random.uniform(low=snrmin, high=snrmax, size=nnoise)
		for qq in range(0,nnoise):
			
			
			# Corrupt signals with noise
			snr = snr_rndvals[qq]
			sig_noisy = np.sqrt( ( sig + (1.0/float(snr))*np.random.randn(sig.size) )**2  +  ( (1.0/float(snr))*np.random.randn(sig.size) )**2 )			
			
			# Fit signals with custom-written code
			pars = vfit(sig_noisy,bvals)
			adc = pars[1]
			akc = pars[2]
			adcdata[nn,qq] = adc
			akcdata[nn,qq] = akc



	## Calculate reference value for infinite SNR
	h = open(sigfile,'rb')
	data = pk.load(h)
	h.close()
	data = data[0]
	sig = data[0]
	bvals = np.loadtxt(bvalfile)
	pars = vfit(sig,bvals)
	adcref = pars[1]
	akcref = pars[2]


	### Plot
	plt.subplot(len(bmax_array),2,plt_count)
	plt.boxplot(np.transpose(adcdata))
	plt.plot([1.0,float(nsub)],[adcref,adcref],'--',linewidth=LW,label='Ref.: 19 b-values, SNR → ∞ ')
	plt.ylabel('$D$  [$\mu$m$^2$/ms]',fontsize=15)
	plt.ylim([0.04,0.25])
	plt.yticks(ticks=[0.05,0.10,0.15,0.20,0.25],fontsize=15)
	plt.xticks(ticks=[1,2,3,4,5],labels=myxlabels,rotation=45,fontsize=13)
	plt.title('{} Max b = {} s/mm$^2$'.format(plt_count_title[plt_count-1],bmax),fontsize=18,fontweight='bold')
	plt.legend(fontsize=15)
	plt.grid()
	plt_count = plt_count + 1

	plt.subplot(len(bmax_array),2,plt_count)
	plt.boxplot(np.transpose(akcdata))
	plt.plot([1.0,float(nsub)],[akcref,akcref],'--',linewidth=LW,label='Ref.: 19 b-values, SNR → ∞ ')
	plt.ylabel('$K$',fontsize=15)
	plt.ylim([-6.0,11.0])
	plt.yticks(ticks=[-4,-2,0,2,4,6,8,10],fontsize=15)
	plt.xticks(ticks=[1,2,3,4,5],labels=myxlabels,rotation=45,fontsize=13)
	plt.title('{} Max b = {} s/mm$^2$'.format(plt_count_title[plt_count-1],bmax),fontsize=18,fontweight='bold')
	plt.legend(fontsize=15)
	plt.grid()
	plt_count = plt_count + 1


# Set up figure
plt.subplots_adjust(left=0.10, bottom=0.11, right=0.898, top=0.90, wspace=0.257, hspace=0.533)
fig = plt.gcf()
fig.set_size_inches([15.0,12.0])

# Show
plt.savefig('script12b_plot.600dpi.png',dpi=600)
plt.show()
print('    ... done ')
print('')

			
