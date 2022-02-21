# Get patch-wise histological parametric maps for each of two MRI slices of two specimens
filelist=("wt/histology/slice1.txt" "wt/histology/slice2.txt" "pdx/histology/slice1.txt" "pdx/histology/slice2.txt")

for FILE in ${filelist[@]}
do

	
	infile=$FILE".csv"
	outroot=$FILE
	python getPatchMapFromQuPath.py $infile 25000 25000 348.8 348.8 1000.0 $outroot --szmax 42.0
	


done

