### Calculate cell size and diffusivity from ADC and Kurtosis (PolyMap method)

# Load packages
import numpy as np
import nibabel as nib
import pickle as pk
import AdcAkcFit as dfit
import time

# Input information
mriroot='.'     # Root directory

scandir=[]     # List of scanned samples
scandir.append('wt/NIFTI_preproc')
scandir.append('pdx/NIFTI_preproc')
snr_values=np.array([68,100])   # SNR values per sample
sigfile='sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift.nii'          # MRI signals
bvalfile='sl_TE45D30_bmin1700.bval'                                             # b-values
outdir='sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC'           # Output folder



## Process all scans
for ss in range(0,len(scandir)):


	# Print information
	tic = time.time()
	print('')
	print('++++++++++++++++++++++++++++++++++++++++++')
	print('         {}'.format(scandir[ss]))
	print('++++++++++++++++++++++++++++++++++++++++++')
	print('')
	
	# Load SNR-dependent estimators
	outstr='polymapSNR{}'.format(snr_values[ss])                                           # Root name for output files
	Lest='exvivomri_d10.0_D30.0/LfitPoly_d10.0_D30.0.sig.icell.SNR{}.bmin1700.bin'.format(snr_values[ss])      # L(ADC,AKC) estimator
	D0est='exvivomri_d10.0_D30.0/D0fitPoly_d10.0_D30.0.sig.icell.SNR{}.bmin1700.bin'.format(snr_values[ss])    # D0(ADC,AKC) estimator
	h = open(Lest,'rb')
	Lfit = pk.load(h)
	h.close()
	h = open(D0est,'rb')
	D0fit = pk.load(h)
	h.close()


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
	
	# Allocate output maps
	Adcmap = np.zeros((mri_size[0],mri_size[1],mri_size[2]))
	Akcmap = np.zeros((mri_size[0],mri_size[1],mri_size[2]))
	Lmap = np.zeros((mri_size[0],mri_size[1],mri_size[2]))
	D0map = np.zeros((mri_size[0],mri_size[1],mri_size[2]))
	S0map = np.zeros((mri_size[0],mri_size[1],mri_size[2]))
	SigPredMap = np.zeros((mri_size[0],mri_size[1],mri_size[2],mri_size[3]))
	
	# Fit (ADC,AKC) and map those values to (D0,L)
	for ii in range(0,mri_size[0]):
		for jj in range(0,mri_size[1]):
			for kk in range(0,mri_size[2]):
			
				if(mask_data[ii,jj,kk]):
			
			
					# Get (ADC,AKC)
					mysig = sig_data[ii,jj,kk,:]
					outpars = dfit.akcfit(mysig, bvals)
					s0 = outpars[0]
					D = outpars[1]
					K = outpars[2]
						
					# Map ADC and AKC to cell size L and cell diffusivity D0
					Lval = Lfit[0] + Lfit[1]*D + Lfit[2]*K + Lfit[3]*D*D + Lfit[4]*K*K + Lfit[5]*D*K + Lfit[6]*D*D*K + Lfit[7]*D*K*K + Lfit[8]*D*D*D +  Lfit[9]*K*K*K
					D0val = D0fit[0] + D0fit[1]*D + D0fit[2]*K + D0fit[3]*D*D + D0fit[4]*K*K + D0fit[5]*D*K + D0fit[6]*D*D*K + D0fit[7]*D*K*K + D0fit[8]*D*D*D +  D0fit[9]*K*K*K

					if(Lval<=0):
						Lval = np.nan
					if(D0val<=0):
						D0val = np.nan
					
					# Store results
					Adcmap[ii,jj,kk] = D
					Akcmap[ii,jj,kk] = K
					Lmap[ii,jj,kk] = Lval
					D0map[ii,jj,kk] = D0val
					S0map[ii,jj,kk] = s0
					bvrescale = bvals/1000.0
					SigPredMap[ii,jj,kk,:] = s0*np.exp( -bvrescale*D + (1.0/6.0)*D*D*K*bvrescale*bvrescale )

	# Save output NIFTI: ADC
	ADCpath = '{}/{}/{}/{}_{}'.format(mriroot,scandir[ss],outdir,outstr,'ADCum2ms.nii')
	ADC_obj = nib.Nifti1Image(Adcmap,buffer_affine,buffer_header)
	nib.save(ADC_obj, ADCpath)

	# Save output NIFTI: AKC
	AKCpath = '{}/{}/{}/{}_{}'.format(mriroot,scandir[ss],outdir,outstr,'AKCum2ms.nii')
	AKC_obj = nib.Nifti1Image(Akcmap,buffer_affine,buffer_header)
	nib.save(AKC_obj, AKCpath)
				
	# Save output NIFTI: cell size L
	Lpath = '{}/{}/{}/{}_{}'.format(mriroot,scandir[ss],outdir,outstr,'Lum.nii')
	L_obj = nib.Nifti1Image(Lmap,buffer_affine,buffer_header)
	nib.save(L_obj, Lpath)		
	
	# Save output NIFTI: cell diffusivity D0
	D0path = '{}/{}/{}/{}_{}'.format(mriroot,scandir[ss],outdir,outstr,'D0um2ms.nii')
	D0_obj = nib.Nifti1Image(D0map,buffer_affine,buffer_header)
	nib.save(D0_obj, D0path)
	
	# Save output NIFTI: S0
	S0path = '{}/{}/{}/{}_{}'.format(mriroot,scandir[ss],outdir,outstr,'S0.nii')
	S0_obj = nib.Nifti1Image(S0map,buffer_affine,buffer_header)
	nib.save(S0_obj, S0path)
	
	# Save output NIFTI: signal prediction
	SigEstpath = '{}/{}/{}/{}_{}'.format(mriroot,scandir[ss],outdir,outstr,'SigPred.nii')
	SigEst_obj = nib.Nifti1Image(SigPredMap,buffer_affine,buffer_header)
	nib.save(SigEst_obj, SigEstpath)
	
	# Print how long it took
	toc = time.time()
	print('     ------> it took {} seconds'.format(toc-tic))
	print('')

