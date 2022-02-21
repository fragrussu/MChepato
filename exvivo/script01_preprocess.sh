#### Run denoising and Gibbs ringing mitigation on the PDX case

DATADIR="."
FSLOUTPUTTYPE=NIFTI

for DIR in $DATADIR"/pdx/NIFTI_preproc" $DATADIR"/wt/NIFTI_preproc"
do

	python $MRIPYTOOLS"/runMPPCA.py" $DIR"/sl_TE31D15_TE45D15_TE65D15_TE45D30_TE65D30_TR2250_TR3300.nii" $DIR"/sl_TE31D15_TE45D15_TE65D15_TE45D30_TE65D30_TR2250_TR3300_mppca7x7x3" --kernel 7,7,3 --nthread 4 --nsa 1   # Denoise with MP-PCA and mitigate noise floor
	fslmaths $DIR"/sl_TE31D15_TE45D15_TE65D15_TE45D30_TE65D30_TR2250_TR3300_mppca7x7x3_den_unbias.nii" -nan $DIR"/sl_TE31D15_TE45D15_TE65D15_TE45D30_TE65D30_TR2250_TR3300_mppca7x7x3_den_unbias.nii"             # Remove NaNs
	python $MRIPYTOOLS"/unring.py" $DIR"/sl_TE31D15_TE45D15_TE65D15_TE45D30_TE65D30_TR2250_TR3300_mppca7x7x3_den_unbias.nii" $DIR"/sl_TE31D15_TE45D15_TE65D15_TE45D30_TE65D30_TR2250_TR3300_mppca7x7x3_den_unbias_unring.nii" --np 18 --axis 2   # Gibbs unringing



done
