## Plot histological and MRI parametric maps

import sys
import matplotlib 
from matplotlib import pyplot as plt
import matplotlib as mpl
import nibabel as nib
import numpy as np
from skimage import morphology as mph 

# Add plotting tools
sys.path.insert(0,'/home/fgrussu/lib/python/MRItools/tools')  # path to MRItools
import plottools as ptls



## Input information
rootdir="."
specdir=["pdx","wt"]
snrval=[100,68]
zslice=[1,2]
xmin = [10,10] 
xmax = [80,80] 
ymin = [10,25] 
ymax = [100,115]

# Colourmap
mycmapname = 'viridis'
cmap = mpl.cm.get_cmap(mycmapname)
mycolours = cmap.colors
mymap = np.zeros((len(mycolours),3))
for qq in range(0,len(mycolours)):
	mymap[qq,0] = mycolours[qq][0]
	mymap[qq,1] = mycolours[qq][1]
	mymap[qq,2] = mycolours[qq][2]


# Metric range
Lpoly_min = 5.0
Lpoly_max = 55.0
Lsig_min = 13.0
Lsig_max = 30.0
Lhisto_min = 13.0
Lhisto_max = 26.0
D0poly_min = 0.05
D0poly_max = 1.95
D0sig_min = 0.05
D0sig_max = 1.95

## Load MRI and histology maps for all specimens
img_all = np.array([])
for ss in range(0,len(specdir)):


	# Load anatomical MRI information
	buffnii = nib.load('{}/{}/NIFTI_preproc/sl_TE31D15_TE45D15_TE65D15_TE45D30_TE65D30_mppca7x7x3_den_unbias_unring_undrift.nii'.format(rootdir,specdir[ss]))
	buffnii = buffnii.get_fdata()
	b0 = buffnii[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss],0]
	b0 = np.fliplr(np.rot90(b0))
	del buffnii
	
	buffnii = nib.load('{}/{}/NIFTI_preproc/sl_histo2mri_ODeosin_mask.nii'.format(rootdir,specdir[ss]))
	buffnii = buffnii.get_fdata()
	mask = buffnii[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	mask = np.fliplr(np.rot90(mask))
	mask = mph.binary_erosion(mask, selem=None, out=None)
	mask_empty = 0.0*mask
	del buffnii
	
	# Load histology information
	buffnii = nib.load('{}/{}/NIFTI_preproc/sl_histo2mri_Lum.nii'.format(rootdir,specdir[ss]))
	buffnii = buffnii.get_fdata()
	Lhisto = buffnii[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	Lhisto = np.fliplr(np.rot90(Lhisto))
	del buffnii
	
	# Load MRI parametric maps
	buffnii = nib.load('{}/{}/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR{}_Lum.nii'.format(rootdir,specdir[ss],snrval[ss]))
	buffnii = buffnii.get_fdata()
	Lpoly = buffnii[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	Lpoly = np.fliplr(np.rot90(Lpoly))
	
	buffnii = nib.load('{}/{}/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR{}_D0um2ms.nii'.format(rootdir,specdir[ss],snrval[ss]))
	buffnii = buffnii.get_fdata()
	D0poly = buffnii[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	D0poly = np.fliplr(np.rot90(D0poly))

	buffnii = nib.load('{}/{}/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_Lum.nii'.format(rootdir,specdir[ss]))
	buffnii = buffnii.get_fdata()
	Lsigfit = buffnii[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	Lsigfit = np.fliplr(np.rot90(Lsigfit))
	
	buffnii = nib.load('{}/{}/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_D0um2ms.nii'.format(rootdir,specdir[ss]))
	buffnii = buffnii.get_fdata()
	D0sigfit = buffnii[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	D0sigfit = np.fliplr(np.rot90(D0sigfit))


	# Overlay
	b0_rgb = ptls.overlay_qmri_over_anatomical(b0,mask_empty,Lpoly,Lpoly_min,Lpoly_max,mymap)
	Lhisto_rgb = ptls.overlay_qmri_over_anatomical(b0,mask,Lhisto,Lhisto_min,Lhisto_max,mymap)
	Lpoly_rgb = ptls.overlay_qmri_over_anatomical(b0,mask,Lpoly,Lpoly_min,Lpoly_max,mymap)
	Lsigfit_rgb = ptls.overlay_qmri_over_anatomical(b0,mask,Lsigfit,Lsig_min,Lsig_max,mymap)
	D0poly_rgb = ptls.overlay_qmri_over_anatomical(b0,mask,D0poly,D0poly_min,D0poly_max,mymap)
	D0sigfit_rgb = ptls.overlay_qmri_over_anatomical(b0,mask,D0sigfit,D0sig_min,D0sig_max,mymap)
	
	
	# Stack data
	img_subj = np.concatenate((b0_rgb,Lhisto_rgb,Lpoly_rgb,Lsigfit_rgb,D0poly_rgb,D0sigfit_rgb),axis=1)
	
	try:
		img_all = np.concatenate((img_all,img_subj),axis=0)
	
	except:
		img_all = np.copy(img_subj)
	

## Add black space
img_all = np.concatenate( (  np.uint8(np.zeros((30,img_all.shape[1],3))) , img_all, np.uint8(np.zeros((30,img_all.shape[1],3)))  ) , axis=0)
img_all = np.concatenate( (  np.uint8(np.zeros((img_all.shape[0],40,3))) , img_all  ) , axis=1)



## Show
fig, ax = plt.subplots()
shw = ax.imshow(img_all)
plt.xticks([])
plt.yticks([])

## Add a colorbar
bar1 = plt.colorbar(mappable=mpl.cm.ScalarMappable(cmap=mycmapname),ax=ax,orientation='vertical',shrink=0.5,aspect=5)

# Make background black
fig.set_facecolor((0.0, 0.0, 0.0))
fig.set_size_inches(20,20)
plt.savefig('./script11_ShowMaps_600dpi.png',dpi=600)
plt.show()



