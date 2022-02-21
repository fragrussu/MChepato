PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      echo ""
      echo "Usage: getCSVfromQuPath.sh INFILE OUTFILE"
      echo ""
      echo "This program re-formats a cell detection .txt file provided by the QuPath"
      echo "histology analysis software to a CSV-formatted file."
      echo ""
      echo "* INFILE: path of the .txt file as provided by QuPath cell detection"
      echo "* OUTFILE: path of the output CSV file (it contains a header with variable names)"
      echo ""
      echo ""
      echo "OPTIONS"
      echo ""
      echo "  -h, --help                 print this help manual"
      echo ""
      echo "Author:" 
      echo "- Francesco Grussu <fgrussu@vhio.net>"
      echo ""
      exit	
      shift
      ;;
    *) 
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done

# Get positional arguments
eval set -- "$PARAMS"
infile=$1
outfile=$2

# Print feedback
echo "************************************************"
echo "             getCSVfromQuPath.sh                "
echo "************************************************"
echo ""
echo "* Input file: "$infile 
echo "* Output file: "$outfile
echo ""

# Move data to RAM
bufferfolder="/dev/shm/__"$RANDOM"_"$RANDOM    # Temporary stuff will be stored in RAM
bufferin=$bufferfolder"/in.txt"
bufferout=$bufferfolder"/out.csv"
mkdir $bufferfolder
cp $infile $bufferin


# Loop through lines (discard the first line, which is a header)
echo "   ... counting lines"
nlines=`cat $bufferin | wc -l`
echo "                      --> it has "$nlines" lines"
echo ""
echo "   ... skipping first line (header)"
headerstr="Centroid_X_um,Centroid_Y_um,Nucleus_Area,Nucleus_Perimeter,Nucleus_Circularity,Nucleus_Max_caliper,Nucleus_Min_caliper,Nucleus_Eccentricity,Nucleus_Hematoxylin_OD_mean,Nucleus_Hematoxylin_OD_sum,Nucleus_Hematoxylin_OD_std,Nucleus_Hematoxylin_OD_max,Nucleus_Hematoxylin_OD_min,Nucleus_Hematoxylin_OD_range,Nucleus_Eosin_OD_mean,Nucleus_Eosin_OD_sum,Nucleus_Eosin_OD_std,Nucleus_Eosin_OD_max,Nucleus_Eosin_OD_min,Nucleus_Eosin_OD_range,Cell_Area,Cell_Perimeter,Cell_Circularity,Cell_Max_caliper,Cell_Min_caliper,Cell_Eccentricity,Cell_Hematoxylin_OD_mean,Cell_Hematoxylin_OD_std,Cell_Hematoxylin_OD_max,Cell_Hematoxylin_OD_min,Cell_Eosin_OD_mean,Cell_Eosin_OD_std,Cell_Eosin_OD_max,Cell_Eosin_OD_min,Cytoplasm_Hematoxylin_OD_mean,Cytoplasm_Hematoxylin_OD_std,Cytoplasm_Hematoxylin_OD_max,Cytoplasm_Hematoxylin_OD_min,Cytoplasm_Eosin_OD_mean,Cytoplasm_Eosin_OD_std,Cytoplasm_Eosin_OD_max,Cytoplasm_Eosin_OD_min,NucleusCell_area_ratio"
echo $headerstr > $bufferout
echo ""
echo "   ... processing additional lines"
echo ""
for i in `seq 2 $nlines`
do

	echo "         - line "$i"/"$nlines
	wordcount=1
		
	# Create a CSV-like string for the current line
	linestr=""
	for WORD in `cat $bufferin | head -n $i | tail -n 1`
	do
		
		# Discard the first 4 words
		if [ $wordcount -le 4 ]
		then
			donothing=1
			
		# Add all other words one after another, using a comma to separate them
		else
			linestr=$linestr","$WORD
		fi
		let wordcount=wordcount+1
	done	
		
	# Append to the output file
	len_linestr=${#linestr} # String length
	linestr=`echo $linestr | cut -c2-$len_linestr` # Remove the first character, as it is an unwanted comma
	echo $linestr >> $bufferout
	echo ""
		
done

echo ""
echo "   ... cleaning up"
echo ""
cp $bufferout $outfile
rm -r -f $bufferfolder

echo ""
echo "   ... done"
echo ""
	



