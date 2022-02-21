import numpy as np
import nibabel as nib
import sys
from skimage import morphology as mph
from matplotlib import pyplot as plt
import matplotlib as mpl

# Add plotting tools
sys.path.insert(0,'<PATH TO MRI TOOLS FROM https://github.com/fragrussu/MRItools>')  
import plottools as ptls


zslice = [2,1]
xval = [45,57]
yval = [84,84] 
mytitles = ['C) Wild type','D) PDX']

signals = []
signals.append('wt/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift.nii')
signals.append('pdx/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift.nii')

bvals = []
bvals.append('wt/NIFTI_preproc/sl_TE45D30_bmin1700.bval')
bvals.append('pdx/NIFTI_preproc/sl_TE45D30_bmin1700.bval')

adcpolymap = []
adcpolymap.append('wt/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR68_ADCum2ms.nii')
adcpolymap.append('pdx/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR100_ADCum2ms.nii')

akcpolymap = []
akcpolymap.append('wt/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR68_AKCum2ms.nii')
akcpolymap.append('pdx/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR100_AKCum2ms.nii')

s0polymap = []
s0polymap.append('wt/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR68_S0.nii')
s0polymap.append('pdx/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR100_S0.nii')

predpolymap = []
predpolymap.append('wt/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR68_SigPred.nii')
predpolymap.append('pdx/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR100_SigPred.nii')

lsigfit = []
lsigfit.append('wt/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_Lum.nii')
lsigfit.append('pdx/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_Lum.nii')

d0sigfit = []
d0sigfit.append('wt/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_D0um2ms.nii')
d0sigfit.append('pdx/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_D0um2ms.nii')

s0sigfit = []
s0sigfit.append('wt/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_S0.nii')
s0sigfit.append('pdx/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_S0.nii')

predsigfit = []
predsigfit.append('wt/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_SigPred.nii')
predsigfit.append('pdx/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_SigPred.nii')
Nsamples = len(signals)

## Loop through samples for plotting
for ss in range(0,Nsamples):

	### Load signals
	sig = nib.load(signals[ss])
	sig = sig.get_fdata()

	# Load b-values
	bv = np.loadtxt(bvals[ss])/1000.0    # Get b-values in ms/um2
	delta = 10.0
	Delta = 30.0
	
	# Load metrics from PolyMap and SigFit
	adcpoly = nib.load(adcpolymap[ss])
	adcpoly = adcpoly.get_fdata()
	
	akcpoly = nib.load(akcpolymap[ss])
	akcpoly = akcpoly.get_fdata()

	s0poly = nib.load(s0polymap[ss])
	s0poly = s0poly.get_fdata()
	
	lsig = nib.load(lsigfit[ss])
	lsig = lsig.get_fdata()

	d0sig = nib.load(d0sigfit[ss])
	d0sig = d0sig.get_fdata()
	
	s0sig = nib.load(s0sigfit[ss])
	s0sig = s0sig.get_fdata()

	
	## Get signal and metrics within voxel
	sigarray = sig[xval[ss],yval[ss],zslice[ss]]
	ADCp = adcpoly[xval[ss],yval[ss],zslice[ss]]
	AKCp = akcpoly[xval[ss],yval[ss],zslice[ss]]
	S0p = s0poly[xval[ss],yval[ss],zslice[ss]]
	Ls = lsig[xval[ss],yval[ss],zslice[ss]]
	D0s = d0sig[xval[ss],yval[ss],zslice[ss]]
	S0s = s0sig[xval[ss],yval[ss],zslice[ss]]
	a0 = 0.0013416743239695044   # Constants to compute ADC
	a1 = 1.2593797077696043e-05  # Constants to compute ADC
	ADCs = a0*(Ls**4)/(D0s*delta*( Delta - delta/3.0 )) - a1*(Ls**6)/(D0s*D0s*delta*delta*( Delta - delta/3.0 ))

	## Synthesise signals
	bvdense = np.linspace(np.min(bv),np.max(bv),100)
	predpolymap = S0p*np.exp(-bvdense*ADCp + (1.0/6.0)*ADCp*ADCp*AKCp*bvdense*bvdense)
	predsigfit = S0s*np.exp(-bvdense*ADCs)
	
	## Plot
	plt.subplot(1,2,ss+1)
	plt.plot(bvdense*1000,np.log(predpolymap), linewidth=3, label='$PolyMap$ fitting')
	plt.plot(bvdense*1000,np.log(predsigfit), linewidth=3, label='$SigFit$ fitting')
	plt.plot(bv*1000,np.log(sigarray),'ko', markersize=10, label='Measurements (log scale)')
	plt.xticks(bv*1000, fontsize=13)
	plt.yticks(fontsize=13)
	plt.ylabel('$ln(s)$', fontsize=15)
	plt.xlabel('$b$   [s/mm$^2$]', fontsize=15)
	plt.title(mytitles[ss], fontsize=18, fontweight='bold')
	plt.legend(fontsize=13)
	

fig = plt.gcf()
fig.set_size_inches([11.0,5.0])
plt.subplots_adjust(left=0.093, bottom=0.145, right=0.971, top=0.88, wspace=0.32, hspace=0.2)
plt.savefig('script13_signalfitting_600dpi.png',dpi=600)
plt.show()



