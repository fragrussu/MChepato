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
Lmin = 13.0
Lmax = 35.0

mycmapname = 'viridis'
d0list = np.array([0.5, 0.75, 1.0, 1.25, 1.5])
Ndiff = d0list.size

scanlist = []
scanlist.append('wt/NIFTI_preproc')
scanlist.append('pdx/NIFTI_preproc')
Nsamples = len(scanlist) 

for ss in range(0,Nsamples):

	### Pint sample info
	print('****  Sample {}'.format(scanlist[ss]))

	### Load data
	# Mask
	mk = nib.load('{}/sl_histo2mri_ODeosin_mask.nii'.format(scanlist[ss]))
	mk = mk.get_fdata()
	mk = mk[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	mk = np.fliplr(np.rot90(mk))
	mk = mph.binary_erosion(mk, selem=None, out=None)
	
	# Reference L
	lsig = nib.load('{}/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_Lum.nii'.format(scanlist[ss]))
	lsig = lsig.get_fdata()
	lsig = lsig[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
	lsig = np.fliplr(np.rot90(lsig))
	lsig[mk==0] = np.nan       # Store NaN on the background for visualisation purposes
	
	# Loop through different fixed diffusivity values
	for qq in range(0,Ndiff):
	
		# Load L estimated while fixing D0
		lsigfix = nib.load('{}/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_Lum_D0fixed{}um2ms.nii'.format(scanlist[ss],d0list[qq]))
		lsigfix = lsigfix.get_fdata()
		lsigfix = lsigfix[xmin[ss]:xmax[ss],ymin[ss]:ymax[ss],zslice[ss]]
		lsigfix = np.fliplr(np.rot90(lsigfix))
		lsigfix[mk==0] = np.nan   # Store NaN on the background for visualisation purposes
		
		# Concatenate maps obtained fixing D0 on the same sample
		if(qq==0):
			lsigfix_frame = np.copy(lsigfix)
		else:
			lsigfix_frame = np.concatenate((lsigfix_frame,lsigfix),axis=1)
			
	# Concatenate all maps from the same sample
	myout = np.concatenate((lsig,lsigfix_frame),axis=1)
	
	# Concatenate samples
	if(ss==0):
		allmaps = np.copy(myout)
	else:
		allmaps = np.concatenate((allmaps,myout),axis=0)


							
# Show
fig, ax = plt.subplots()
shw = ax.imshow(allmaps,interpolation='None',vmin=Lmin,vmax=Lmax)
current_cmap = mpl.cm.get_cmap()     # Use black for NaNs (i.e. for background)
current_cmap.set_bad(color='black')  # Use black for NaNs (i.e. for background)
plt.xticks([])
plt.yticks([])
fig.set_facecolor((0.0, 0.0, 0.0))
fig.set_size_inches(20,20)
bar1 = plt.colorbar(mappable=mpl.cm.ScalarMappable(cmap=mycmapname),ax=ax,orientation='vertical',shrink=0.5,aspect=5)
plt.savefig('./script14_SigFitD0fixed_Lmin{:.2f}_Lmax{:.2f}_600dpi.png'.format(Lmin,Lmax),dpi=600)
plt.show()
	


