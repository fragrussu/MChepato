### Calculate SNR after denoising and estimate mean/std of SNR within the sample

## I/O information
export FSLOUTPUTTYPE=NIFTI
rootdir="."

## SNR calculation
for sample in "wt" "pdx"
do
	
	# Calculate SNR
	liver=$rootdir"/"$sample"/NIFTI_preproc/sl_livermask.nii"
	dwi=$rootdir"/"$sample"/NIFTI_preproc/sl_TE31D15_TE45D15_TE65D15_TE45D30_TE65D30_TR2250_TR3300_mppca7x7x3_den.nii"
	sigma=$rootdir"/"$sample"/NIFTI_preproc/sl_TE31D15_TE45D15_TE65D15_TE45D30_TE65D30_TR2250_TR3300_mppca7x7x3_sigma.nii"
	snr=$rootdir"/"$sample"/NIFTI_preproc/sl_TE31D15_TE45D15_TE65D15_TE45D30_TE65D30_TR2250_TR3300_mppca7x7x3_snrmap.nii"
	fslmaths $dwi -div $sigma $snr

	# Print information
	echo "+++++++++++++++++++++++++++++++++++++++++++++"
	echo "            Sample "$sample
	echo "+++++++++++++++++++++++++++++++++++++++++++++"
	echo ""
	
	# Evaluate SNR at b = 0 for TE = 45 ms
	fslroi $snr $rootdir"/"$sample"/NIFTI_preproc/dummyB.nii" 10 1
	fslroi $snr $rootdir"/"$sample"/NIFTI_preproc/dummyC.nii" 30 1
	fslmerge -t $rootdir"/"$sample"/NIFTI_preproc/dummyD.nii" $rootdir"/"$sample"/NIFTI_preproc/dummyB.nii" $rootdir"/"$sample"/NIFTI_preproc/dummyC.nii"
	fslmaths $rootdir"/"$sample"/NIFTI_preproc/dummyD.nii" -Tmean $rootdir"/"$sample"/NIFTI_preproc/dummyD.nii"
	snrmean45=`fslstats $rootdir"/"$sample"/NIFTI_preproc/dummyD.nii" -k $liver -n -m`
	snrmed45=`fslstats $rootdir"/"$sample"/NIFTI_preproc/dummyD.nii" -k $liver -n -p 50`
	echo "    * SNR at TE = 45ms --> mean "$snrmean45", median "$snrmed45
	
	# Clean it up
	rm $rootdir"/"$sample"/NIFTI_preproc/dummyA.nii" $rootdir"/"$sample"/NIFTI_preproc/dummyB.nii" $rootdir"/"$sample"/NIFTI_preproc/dummyC.nii" $rootdir"/"$sample"/NIFTI_preproc/dummyD.nii"
	echo ""
	
done


