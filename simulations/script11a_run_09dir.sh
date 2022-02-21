# Study the effect of the number of gradient directions: use 9 directions per b-value

mkdir "results_script11"
ndir="9"
python script11_synDKI_OneD0andLVal_SNRinf.py 20.0 75.0 100.0 2000.0 7 "dir"$ndir".txt" "results_script11/dir"$ndir
