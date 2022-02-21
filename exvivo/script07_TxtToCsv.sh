## Convert QuPath detection measurement files (.txt) to CSV (.csv)

filelist=("wt/histology/slice1" "wt/histology/slice2" "pdx/histology/slice1" "pdx/histology/slice2")

for FILE in ${filelist[@]}
do

	# Get file
	infile=$FILE".txt"
	outfile=$FILE".csv"
	
	./getCSVfromQuPath.sh $infile $outfile
	


done





