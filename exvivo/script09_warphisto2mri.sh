### Coregister  manual outlines of tissue in both histology and MRI, and then apply the warp to histological maps

export FSLOUTPUTTYPE=NIFTI
filelist=("wt" "pdx")

for FILE in ${filelist[@]}
do

	# Stuff common to both histological slices
	histodir="./"$FILE
	mridir="./"$FILE"/NIFTI_preproc"
	mrifile=$mridir"/sl_livermask.nii"

	# First slice
	drivesegmri=$mridir"/sl_livermask_drivecoreg.nii"
	histofile=$histodir"/slice1_livermask.nii"
	histomaps=$histodir"/slice1_ODeosin.nii,"$histodir"/slice1_Lum.nii,"$histodir"/slice1_Cellsmm2.nii,"$histodir"/slice1_FCellPatch.nii"
	outroot=$histodir"/slice1"
	merger1=$outroot   # To merge all slices to one NIFTI
	zslice="1"

	python warpPatchHisto2MRI.py $drivesegmri $zslice $histofile $outroot --histo_maps $histomaps
	mv -v $outroot"_histo2mri_map1.nii" $outroot"_histo2mri_ODeosin.nii"
	mv -v $outroot"_histo2mri_map2.nii" $outroot"_histo2mri_Lum.nii"
	mv -v $outroot"_histo2mri_map3.nii" $outroot"_histo2mri_Cellsmm2.nii"
	mv -v $outroot"_histo2mri_map4.nii" $outroot"_histo2mri_FCellPatch.nii"


	# Second slice
	drivesegmri=$mridir"/sl_livermask_drivecoreg.nii"
	histofile=$histodir"/slice2_livermask.nii"
	histomaps=$histodir"/slice2_ODeosin.nii,"$histodir"/slice2_Lum.nii,"$histodir"/slice2_Cellsmm2.nii,"$histodir"/slice2_FCellPatch.nii"
	outroot=$histodir"/slice2"
	merger2=$outroot   # To merge all slices to one NIFTI
	zslice="2"

	python warpPatchHisto2MRI.py $drivesegmri $zslice $histofile $outroot --histo_maps $histomaps
	mv -v $outroot"_histo2mri_map1.nii" $outroot"_histo2mri_ODeosin.nii"
	mv -v $outroot"_histo2mri_map2.nii" $outroot"_histo2mri_Lum.nii"
	mv -v $outroot"_histo2mri_map3.nii" $outroot"_histo2mri_Cellsmm2.nii"
	mv -v $outroot"_histo2mri_map4.nii" $outroot"_histo2mri_FCellPatch.nii"

	# Merge the two slices to one 3D NIFTI
	echo ""
	echo "... copying all histological files warped to MRI space to the MRI folder"
	echo ""
	fslmaths $merger1"_histo2mri_ODeosin.nii" -add $merger2"_histo2mri_ODeosin.nii" -mul $mrifile $mridir"/sl_histo2mri_ODeosin.nii"
	fslmaths $merger1"_histo2mri_Lum.nii" -add $merger2"_histo2mri_Lum.nii" -mul $mrifile $mridir"/sl_histo2mri_Lum.nii"
	fslmaths $merger1"_histo2mri_Cellsmm2.nii" -add $merger2"_histo2mri_Cellsmm2.nii" -mul $mrifile $mridir"/sl_histo2mri_Cellsmm2.nii"
	fslmaths $merger1"_histo2mri_FCellPatch.nii" -add $merger2"_histo2mri_FCellPatch.nii" -mul $mrifile $mridir"/sl_histo2mri_FCellPatch.nii"
	fslmaths $merger1"_histo2mri_tissuemask.nii" -add $merger2"_histo2mri_tissuemask.nii" $mridir"/sl_histo2mri_histoutline.nii"


	# Threshold ODeosin to remove areas for which there is tissue in MRI but not in histology
	fslmaths $mridir"/sl_histo2mri_ODeosin.nii" -uthr 0.999 -bin $mridir"/sl_histo2mri_ODeosin_mask.nii"
	
done

