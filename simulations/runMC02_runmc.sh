### Run MCDC diffusion simulator on prism models of hepatocytes for various diffusivities
#

meshdir="../data/perturbed"   # Directory where meshes are stored
spindir="../data/perturbed_randomwalks"    # Directory where random walks will be stored
mkdir $spindir
mkdir "./runMC02_runmc.info"  # Directory where configuration files will be stored

echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "             Random walks with MCDC simulator                 "
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo ""
echo ""

# Loop over sizes
for mysize in 6.5 8.0 9.5 11.0 12.5 14.0 15.5 17.0 18.5 20.0 21.5 23.0 24.5 26.0 27.5 29.0 30.5 32.0 33.5 35.0 36.5 38.0 39.5 41.0 42.5 44.0 45.5 47.0 48.5 50.0 51.5 53.0 54.5 56.0 57.5 60.0
do
	# Loop over diffusivity values
	for dval in 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50 1.55 1.60 1.65 1.70 1.75 1.80 1.85 1.90 1.95 2.00 2.05 2.10 2.15 2.20 2.25 2.30 2.35 2.40
	do


		# Loop over basic shapes (4, 5 and 6 sides)
		for myside in 4 5 6
		do


			# Loop over different cell instantiations (perturbations)
			for nn in `seq 0 4`
			do

				# Path of text files
				confile="./script02_runmc.info/mc_sides"$myside"_diam"$mysize"um_fpert0.1_npert"$nn"_dval"$dval"um2ms.conf"
				confile_linediff="./script02_runmc.info/mc_sides"$myside"_diam"$mysize"um_fpert0.1_npert"$nn"_dval"$dval"um2ms.conf.ld"
				confile_lineprefix="./script02_runmc.info/mc_sides"$myside"_diam"$mysize"um_fpert0.1_npert"$nn"_dval"$dval"um2ms.conf.lp"
				confile_lineobst="./script02_runmc.info/mc_sides"$myside"_diam"$mysize"um_fpert0.1_npert"$nn"_dval"$dval"um2ms.conf.lo"
				meshfile=$meshdir"/sides"$myside"_diam"$mysize"um_fpert0.1_npert"$nn".ply"
				spinroot=$spindir"/mc_sides"$myside"_diam"$mysize"um_fpert0.1_npert"$nn"_dval"$dval"um2ms.conf"
				
				# Create configuration file for MCDC simulator
				echo "diffusivity "$dval"e-6" > $confile_linediff
				echo "exp_prefix "$spinroot > $confile_lineprefix
				echo "ply "$meshfile > $confile_lineobst
				echo ""
				echo "  ... creating configuration file:"
				echo "      "$confile
				cat "./runMC02_runmc.confp1.conf" $confile_linediff $confile_lineprefix "./runMC02_runmc.confp2.conf" $confile_lineobst "./runMC02_runmc.confp3.conf" > $confile

				rm $confile_linediff $confile_lineprefix $confile_lineobst
				
				# Run simulation
				MC-CD_Simulator $confile 
				



			done



		done

	done


done



