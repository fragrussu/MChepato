### Calculate cell size and diffusivity from direct fitting of a biophysical model (SigFit method)

# Load packages
import numpy as np
import nibabel as nib
import pickle as pk
import AdcAkcFit as dfit
import getDict as gt


# Input information
mriroot='.'     # Root directory

scandir=[]     # Scan list
scandir.append('wt/NIFTI_preproc')
scandir.append('pdx/NIFTI_preproc')
sigfile='sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift.nii'          # MRI signals
bvalfile='sl_TE45D30_bmin1700.bval'                                             # b-values
outdir='sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC'           # Output folder

# Dictionary size for fitting
dictsize = 60

## Process all scans
for ss in range(0,len(scandir)):


	# Print information
	print('')
	print('++++++++++++++++++++++++++++++++++++++++++')
	print('         {}'.format(scandir[ss]))
	print('++++++++++++++++++++++++++++++++++++++++++')
	print('')

	# Load signal
	sigpath = '{}/{}/{}'.format(mriroot,scandir[ss],sigfile)
	sig_obj = nib.load(sigpath)
	sig_data = sig_obj.get_fdata()
	mri_size = sig_data.shape
	buffer_header = sig_obj.header    # Reference NIFTI header
	buffer_affine = sig_obj.affine    # Reference NIFTI header scanner-to-image affine transformation
	buffer_header.set_data_dtype('float64')   # Change data type to float64 in reference header

	# Load liver mask
	maskpath = '{}/{}/{}'.format(mriroot,scandir[ss],'sl_livermask.nii')
	mask_obj = nib.load(maskpath)
	mask_data = mask_obj.get_fdata()
	
	# Load b-values
	bvals = np.loadtxt('{}/{}/{}'.format(mriroot,scandir[ss],bvalfile))
	bvals = np.array(bvals)
	delta = 10.0 + np.zeros(bvals.size)
	Delta = 30.0 + np.zeros(bvals.size)
	
	# Allocate output maps
	Lmap = np.zeros((mri_size[0],mri_size[1],mri_size[2]))
	D0map = np.zeros((mri_size[0],mri_size[1],mri_size[2]))
	S0map = np.zeros((mri_size[0],mri_size[1],mri_size[2]))
	SigPredMap = np.zeros((mri_size[0],mri_size[1],mri_size[2],mri_size[3]))
	
	# Fit analytical model of intra-cellular diffusion
	for ii in range(0,mri_size[0]):
		for jj in range(0,mri_size[1]):
			for kk in range(0,mri_size[2]):
			
				if(mask_data[ii,jj,kk]):
			
					# Get MRI signal
					mysig = sig_data[ii,jj,kk,:]
					
					# Estimate S0
					outpars = dfit.adcfit(mysig, bvals)
					s0fit = outpars[0]
					S0map[ii,jj,kk] = s0fit
					
					# Normalise signal
					mysig = mysig/mysig[0]
					
					# Get dictionary of synthetic signals
					synsig,lvals,d0vals = gt.syn(bvals,delta,Delta,nwords=dictsize)
					
					# Normalise dictionary
					for rr in range(0,synsig.shape[0]):
						synsig[rr,:] = synsig[rr,:]/synsig[rr,0]
					
					# Compare signals to dictionary
					mysig_mat = np.tile(mysig,(synsig.shape[0],1))    # Matrix of actual MRI measurements
					mse_array = np.nanmean( (mysig_mat - synsig)**2 , axis=1)
					min_idx = np.argmin(mse_array)
					lfit = lvals[min_idx]
					d0fit = d0vals[min_idx]
					Lmap[ii,jj,kk] = lfit
					D0map[ii,jj,kk] = d0fit
					
					
					# Save predicted signals
					a0 = 0.0013416743239695044   # Constants to compute ADC
					a1 = 1.2593797077696043e-05  # Constants to compute ADC
					adcintra = a0*(lfit**4)/(d0fit*delta*( Delta - delta/3.0 )) - a1*(lfit**6)/(d0fit*d0fit*delta*delta*( Delta - delta/3.0 ))
					SigPredMap[ii,jj,kk,:] = s0fit*np.exp(-(bvals/1000.0)*adcintra)
					
				
				
	# Save output NIFTI: cell size L
	Lpath = '{}/{}/{}/{}'.format(mriroot,scandir[ss],outdir,'sigfit_Lum.nii')
	L_obj = nib.Nifti1Image(Lmap,buffer_affine,buffer_header)
	nib.save(L_obj, Lpath)		
	
	# Save output NIFTI: cell diffusivity D0
	D0path = '{}/{}/{}/{}'.format(mriroot,scandir[ss],outdir,'sigfit_D0um2ms.nii')
	D0_obj = nib.Nifti1Image(D0map,buffer_affine,buffer_header)
	nib.save(D0_obj, D0path)

	# Save output NIFTI: S0
	S0path = '{}/{}/{}/{}'.format(mriroot,scandir[ss],outdir,'sigfit_S0.nii')
	S0_obj = nib.Nifti1Image(S0map,buffer_affine,buffer_header)
	nib.save(S0_obj, S0path)

	# Save output NIFTI: signal prediction
	SigEstpath = '{}/{}/{}/{}'.format(mriroot,scandir[ss],outdir,'sigfit_SigPred.nii')
	SigEst_obj = nib.Nifti1Image(SigPredMap,buffer_affine,buffer_header)
	nib.save(SigEst_obj, SigEstpath)
	

