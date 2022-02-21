import numpy as np
import nibabel as nib
import sys
from skimage import morphology as mph
from matplotlib import pyplot as plt
import matplotlib as mpl

# Add plotting tools
sys.path.insert(0,'<PATH TO MRI TOOLS FROM https://github.com/fragrussu/MRItools>')  
import plottools as ptls


zslice=[2,1]
xmin = [10,10] 
xmax = [80,80] 
ymin = [25,10] 
ymax = [115,100]



masks = []
masks.append('wt/NIFTI_preproc/sl_histo2mri_ODeosin_mask.nii')
masks.append('pdx/NIFTI_preproc/sl_histo2mri_ODeosin_mask.nii')

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
for ss in range(0,Nsamples):

	### Load data
	sig = nib.load(signals[ss])
	sig = sig.get_fdata()
	sig = sig[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss],:]
	sig = np.fliplr(np.rot90(sig))
	
	mk = nib.load(masks[ss])
	mk = mk.get_fdata()
	mk = mk[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	mk = np.fliplr(np.rot90(mk))
	mk = mph.binary_erosion(mk, selem=None, out=None)
	
	bv = np.loadtxt(bvals[ss])
	
	adcpoly = nib.load(adcpolymap[ss])
	adcpoly = adcpoly.get_fdata()
	adcpoly = adcpoly[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	adcpoly = np.fliplr(np.rot90(adcpoly))
	
	akcpoly = nib.load(akcpolymap[ss])
	akcpoly = akcpoly.get_fdata()
	akcpoly = akcpoly[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	akcpoly = np.fliplr(np.rot90(akcpoly))
	
	s0poly = nib.load(s0polymap[ss])
	s0poly = s0poly.get_fdata()
	s0poly = s0poly[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	s0poly = np.fliplr(np.rot90(s0poly))
	
	predpoly = nib.load(predpolymap[ss])
	predpoly = predpoly.get_fdata()
	predpoly = predpoly[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss],:]
	predpoly = np.fliplr(np.rot90(predpoly))
	

	lsig = nib.load(lsigfit[ss])
	lsig = lsig.get_fdata()
	lsig = lsig[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	lsig = np.fliplr(np.rot90(lsig))

	d0sig = nib.load(d0sigfit[ss])
	d0sig = d0sig.get_fdata()
	d0sig = d0sig[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	d0sig = np.fliplr(np.rot90(d0sig))	

	s0sig = nib.load(s0sigfit[ss])
	s0sig = s0sig.get_fdata()
	s0sig = s0sig[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	s0sig = np.fliplr(np.rot90(s0sig))
	
	predsig = nib.load(predsigfit[ss])
	predsig = predsig.get_fdata()
	predsig = predsig[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss],:]
	predsig = np.fliplr(np.rot90(predsig))
	
	# Loop through volumes: concatenate information
	for qq in range(0,sig.shape[2]):
		if(qq==0):
			mysignal = np.copy(sig[:,:,qq]*mk)
			mypoly = np.copy(predpoly[:,:,qq]*mk)
			mysigfit = np.copy(predsig[:,:,qq]*mk)
		else:
			mysignal = np.concatenate((mysignal,sig[:,:,qq]*mk),axis=1)
			mypoly = np.concatenate((mypoly,predpoly[:,:,qq]*mk),axis=1)
			mysigfit = np.concatenate((mysigfit,predsig[:,:,qq]*mk),axis=1)
			
	# Concatenate signal and prediction
	myout = np.concatenate((mysignal,mypoly,mysigfit),axis=0)
				
	# Convert to RGB
	mymin = np.min(myout)
	mymax = np.max(myout)
	myoutrgb = np.uint8( 255*(myout - np.min(myout)) / (np.max(myout) - np.min(myout)) )
	myoutrgb = np.dstack((myoutrgb,myoutrgb,myoutrgb))
			
	# Show
	fig, ax = plt.subplots()
	shw = ax.imshow(myoutrgb,interpolation='None')
	plt.xticks([])
	plt.yticks([])
	fig.set_facecolor((0.0, 0.0, 0.0))
	fig.set_size_inches(20,20)
	bar1 = plt.colorbar(mappable=mpl.cm.ScalarMappable(cmap='Greys'),ax=ax,orientation='vertical',shrink=0.5,aspect=5)
	plt.savefig('./script12_sample{}_sigmin{:.2f}_sigmax{:.2f}_600dpi.png'.format(ss,mymin,mymax),dpi=600)
	plt.show()
	


