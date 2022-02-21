## Print mean and standard deviation of cell size in MRI and histology

rootdir="."
samplelist=("wt" "pdx")

kk=0
for sample in ${samplelist[@]} 
do

	## Print specimen information
	echo "++++++++++++++++++++++++++++++++++++++"
	echo "       Sample "$sample
	echo "++++++++++++++++++++++++++++++++++++++"
	echo ""

	## Get specimen mask and reference SigFit maps
	tissue_mask=$rootdir"/"$sample"/NIFTI_preproc/sl_histo2mri_ODeosin_mask.nii"
	sigfitL=$rootdir"/"$sample"/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_Lum.nii"
	
	## Compute median and IQR for reference SigFit
	sigfitL_med=`fslstats $sigfitL -k $tissue_mask -n -p 50`
	sigfitL_prcmin=`fslstats $sigfitL -k $tissue_mask -n -p 25`
	sigfitL_prcmax=`fslstats $sigfitL -k $tissue_mask -n -p 75`

	# Print reference SigFit information
	echo "SigFit L: med/iqr "$sigfitL_med" ("$sigfitL_prcmin"; "$sigfitL_prcmax") um"
	echo ""
	
	## Loop through fixed diffusivity values
	for d0fixed in "0.5" "0.75" "1.0" "1.25" "1.5"
	do
	
		# Get SigFit map obtained for fixing D0
		sigfitL=$rootdir"/"$sample"/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_Lum_D0fixed"$d0fixed"um2ms.nii"
	
		# Compute median and IQR for SigFit obtained for fixing D0
		sigfitL_med=`fslstats $sigfitL -k $tissue_mask -n -p 50`
		sigfitL_prcmin=`fslstats $sigfitL -k $tissue_mask -n -p 25`
		sigfitL_prcmax=`fslstats $sigfitL -k $tissue_mask -n -p 75`
		
		# Print SigFit obtained while fixing D0
		echo "SigFit L for D0 = "$d0fixed"um2/ms: med/iqr "$sigfitL_med" ("$sigfitL_prcmin"; "$sigfitL_prcmax") um"
		echo ""
	
	done

	# Increment counter
	let kk=kk+1
	echo ""

done


