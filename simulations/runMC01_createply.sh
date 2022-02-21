## Merge vertices with header and faces information to create viable PLY meshes
#


# data folders
reffolder="../data/ref"
pertfolder="../data/perturbed"

echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "    Creating PLY meshes from header, vertices and faces"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo ""

# loop through shapes (4, 5 and 6 sides)
for myside in 4 5 6
do

	echo "    ... "$myside" sides"

	# reference header and faces
	myheader=$reffolder"/sides"$myside"_diam1um.header"
	myfaces=$reffolder"/sides"$myside"_diam1um.faces"

	# loop through hepatocyte sizes
	for mysize in 5.0 6.5 8.0 9.5 11.0 12.5 14.0 15.5 17.0 18.5 20.0 21.5 23.0 24.5 26.0 27.5 29.0 30.5 32.0 33.5 35.0 36.5 38.0 39.5 41.0 42.5 44.0 45.5 47.0 48.5 50.0 51.5 53.0 54.5 56.0 57.5 60.0
	do


		echo "             size "$mysize" um"
		
		# loop through perturbations
		for nn in `seq 0 4`
		do
		
			echo "                      perturbation "$nn

			# vertex file
			myvert=$pertfolder"/sides"$myside"_diam"$mysize"um_fpert0.1_npert"$nn".vertices"	
			
			# output PLY mesh
			myply=$pertfolder"/sides"$myside"_diam"$mysize"um_fpert0.1_npert"$nn".ply"
			cat $myheader $myvert $myfaces > $myply
		


		done
		


	done


done


