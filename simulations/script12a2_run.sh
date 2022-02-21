# Study the effect of the number of the number of b-values: use 19 b-value with maximum b of 1000 s/mm2

mkdir -v "./results_script12"
bmax="1000"
python script12_synMRI_intracell_SubsetD0andL_snrinf.py "20.0" "75.0" "100.0" $bmax".0" "19" "./results_script12/bmax"$bmax
