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
snrmax = 100       # Maximum SNR value at b = 0
snrmin = 20        # Minimum SNR value at b = 0
nnoise = 1000    # Number of independent noise instantations
myseed = 21072018  # Seed for reproducibility

## Kurtosis comparison: full kurtosis tensor fit and powder-average approximation for different numbers of gradient directions
print('')
print('+++++++  Kurtosis estimation')
print('')
print('')
np.random.seed(myseed)

mddata = np.zeros((ndirs,nnoise))     # Array storing Dipy mean diffusivity for different noise instantiations and number of gradient directions
adcdata = np.zeros((ndirs,nnoise))    # Array storing custom ADC for different noise instantiations and number of gradient directions

mkdata = np.zeros((ndirs,nnoise))     # Array storing Dipy mean kurtosis for different noise instantiations and number of gradient directions
akcdata = np.zeros((ndirs,nnoise))    # Array storing custom AKC for different noise instantiations and number of gradient directions

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
	ntot = gx.size
		
	# Create a gradient table for DiPy and add b = 0 (DiPy requires b = 0 images)
	gall = np.zeros((ntot,3))
	gall[:,0] = gx
	gall[:,1] = gy
	gall[:,2] = gz
	gall = np.concatenate((np.zeros((1,3)),gall),axis=0)
	ball = np.concatenate(([0.0],bvals))
	gtaball = gradient_table(ball,bvecs=gall,b0_threshold=1.0)
	sigall = np.concatenate(([1.0],sig))     # Noise-free signal from the full directional acquisition, including b = 0 images
	
	# Loop through random noise instantiations
	snr_rndvals = np.random.uniform(low=snrmin, high=snrmax, size=nnoise)
	for qq in range(0,nnoise):
		
		
		# Corrupt signals with noise
		snr = snr_rndvals[qq]
		sigall_noisy = np.sqrt( ( sigall + (1.0/float(snr))*np.random.randn(sigall.size) )**2  +  ( (1.0/float(snr))*np.random.randn(sigall.size) )**2 )
		
		# Calculate powder average of noisy signals
		buall = np.unique(ball)
		suall_noisy = np.zeros(buall.size)
		for bb in range(0,buall.size):
			suall_noisy[bb] = np.mean(sigall_noisy[ball==buall[bb]])
			
		# Fit DKI with DiPy
		dkimodel = dki.DiffusionKurtosisModel(gtaball)
		dkifit = dkimodel.fit(sigall_noisy)
		md = dkifit.md
		mk = dkifit.mk(-5.0, 10.0)
		mkdata[nn,qq] = mk
		mddata[nn,qq] = md
		
		# Fit powder-average signal with custom-written code
		pars = vfit(suall_noisy,buall)
		adc = pars[1]
		akc = pars[2]
		akcdata[nn,qq] = akc
		adcdata[nn,qq] = adc

		# Print some feedback
		print('    SNR = {}, {}-dir. protocol {}/{}, repeat {}/{}'.format(snr,dirs[nn],nn+1,ndirs,qq+1,nnoise))
		print('')


# Save results	
print('    ... saving results ')
print('')

hfile = open('./script11b_analyse.md.bin'.format(snr),'wb')
pk.dump(mddata,hfile) 
hfile.close()

hfile = open('./script11b_analyse.adc.bin'.format(snr),'wb')
pk.dump(adcdata,hfile) 
hfile.close()

hfile = open('./script11b_analyse.mk.bin'.format(snr),'wb')
pk.dump(mkdata,hfile) 
hfile.close()

hfile = open('./script11b_analyse.akc.bin'.format(snr),'wb')
pk.dump(akcdata,hfile) 
hfile.close()

print('    ... done ')
print('')

