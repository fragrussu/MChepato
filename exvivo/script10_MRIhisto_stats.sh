## Print distribution of MRI and histological metrics in the two samples

rootdir="."
samplelist=("wt" "pdx")
snrlist=("68" "100")

kk=0
for sample in ${samplelist[@]} 
do

	## Get file names
	tissue_mask=$rootdir"/"$sample"/NIFTI_preproc/sl_histo2mri_ODeosin_mask.nii"
	histoL=$rootdir"/"$sample"/NIFTI_preproc/sl_histo2mri_Lum.nii"
	polyL=$rootdir"/"$sample"/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR"${snrlist[$kk]}"_Lum.nii"
	polyD0=$rootdir"/"$sample"/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/polymapSNR"${snrlist[$kk]}"_D0um2ms.nii"
	sigfitL=$rootdir"/"$sample"/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_Lum.nii"
	sigfitD0=$rootdir"/"$sample"/NIFTI_preproc/sl_TE45D30_bmin1700_mppca7x7x3_den_unbias_unring_undrift_AKC/sigfit_D0um2ms.nii"
	
	
	## Compute median and IQR
	histoL_med=`fslstats $histoL -k $tissue_mask -n -p 50`
	histoL_prcmin=`fslstats $histoL -k $tissue_mask -n -p 25`
	histoL_prcmax=`fslstats $histoL -k $tissue_mask -n -p 75`
	
	polyL_med=`fslstats $polyL -k $tissue_mask -n -p 50`
	polyL_prcmin=`fslstats $polyL -k $tissue_mask -n -p 25`
	polyL_prcmax=`fslstats $polyL -k $tissue_mask -n -p 75`
	
	sigfitL_med=`fslstats $sigfitL -k $tissue_mask -n -p 50`
	sigfitL_prcmin=`fslstats $sigfitL -k $tissue_mask -n -p 25`
	sigfitL_prcmax=`fslstats $sigfitL -k $tissue_mask -n -p 75`
	
	polyD0_med=`fslstats $polyD0 -k $tissue_mask -n -p 50`
	polyD0_prcmin=`fslstats $polyD0 -k $tissue_mask -n -p 25`
	polyD0_prcmax=`fslstats $polyD0 -k $tissue_mask -n -p 75`
	
	sigfitD0_med=`fslstats $sigfitD0 -k $tissue_mask -n -p 50`
	sigfitD0_prcmin=`fslstats $sigfitD0 -k $tissue_mask -n -p 25`
	sigfitD0_prcmax=`fslstats $sigfitD0 -k $tissue_mask -n -p 75`
	
	## Print results
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "                          Sample "$sample
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo ""
	echo "Histology L: med/iqr "$histoL_med" ("$histoL_prcmin"; "$histoL_prcmax") um"
	echo "PolyMap L: med/iqr "$polyL_med" ("$polyL_prcmin"; "$polyL_prcmax") um"
	echo "SigFit L: med/iqr "$sigfitL_med" ("$sigfitL_prcmin"; "$sigfitL_prcmax") um"
	echo ""
	echo "PolyMap D0: med/iqr "$polyD0_med" ("$polyD0_prcmin"; "$polyD0_prcmax") um2/ms"
	echo "SigFit D0: med/iqr "$sigfitD0_med" ("$sigfitD0_prcmin"; "$sigfitD0_prcmax") um2/ms"
	echo ""
	echo ""

	# Increment counter
	let kk=kk+1

done


