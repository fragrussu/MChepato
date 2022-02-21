#### Simulate MRI signals for one specific microstructural scenario by averaging 3 gradient directions

import numpy as np
import argparse
import multiprocessing
import sys
import os
import glob as gb
import pickle as pk
import matplotlib.pyplot as plt
import time
import syn_ivim_distr as syn



if __name__ == "__main__":



	### Print help and parse arguments
	parser = argparse.ArgumentParser(description='Synthesise MRI signals')
	parser.add_argument('gDur', help='gradient duration in ms (small Delta)')
	parser.add_argument('gSep', help='gradient separation in ms (big Delta)')
	parser.add_argument('bmin', help='minimum b-value in s/mm2')
	parser.add_argument('bmax', help='maximum b-value in s/mm2')
	parser.add_argument('nbval', help='number of linearly spaced b-values')
	parser.add_argument('outdir', help='output directory')
	args = parser.parse_args()
	

	### Get MRI sequence information
	mybmin = float(args.bmin)
	mybmax =  float(args.bmax)
	myDelta = float(args.gSep)
	mydelta = float(args.gDur)
	Nbval = int(args.nbval)
	outdir = args.outdir
	
	# Create output folder
	try:
		os.mkdir(outdir)
	except:
		do_nothing = True

	# Simulation information
	outstring = '{}/syn_d{}_D{}_bmin{}_bmax{}_Nbval{}'.format(outdir,mydelta,myDelta,mybmin,mybmax,Nbval)     # Output file string
	Nspins = 1000 # number of spins
	Ntime = 3000  # number of time steps
	Nsteps = Ntime + 1  # number of simulation interations
	Tdur = 0.140  # total simulation duration in sec
	trajdir='/nfs/Data/Disk2/clfgrussu/FGrussu/fgrussu/simulations/data/meshes/prism/perturbed_randomwalks'    # Directory where random walks are stored
	cellsize = np.array([11.0,12.5,14.0,15.5,17.0])
	dval = np.array([2.20,2.25,2.30,2.35,2.40])
	Nwindow = 2  # Sliding window width in the discrete cell size / diffusivity grid

	
	# MRI sequence information
	bvalseq = np.round(np.linspace(mybmin,mybmax,Nbval))   # List of b-value in sec/mm2
	deltaseq = mydelta*np.ones(Nbval)     # List of delta in msec
	Deltaseq = myDelta*np.ones(Nbval)     # List of delta in msec


	# Generate list of microstructures
	ustruct_list = []
	Ncells = cellsize.size
	Ndiff = dval.size
	for cs in range(Nwindow,Ncells-Nwindow):
		for dv in range(Nwindow,Ndiff-Nwindow):		
			ustruct_list.append(np.array([cellsize[cs],dval[dv]]))   # Each element of list contains: [cell_size, diffusivity]


	# Save list of microstructures (centres of sliding windows)
	hfile = open('{}.ustruct.bin'.format(outstring),'wb')
	pk.dump(ustruct_list,hfile) 
	hfile.close()

	# Run computing
	processinglist = []   # List of things to process
	print('++++++++ searching files: ')
	print('')
	for ii in range(0,len(ustruct_list)):
		# Measure how long it takes
		tic = time.time()
		# Find files with trajectories (study sliding windows of 7 cell sizes and 7 diffusivities
		Lcentre = ustruct_list[ii][0]
		Lcentreidx = np.asscalar( np.squeeze( np.where(cellsize==Lcentre) ) )
		D0centre = ustruct_list[ii][1]
		D0centreidx = np.asscalar( np.squeeze( np.where(dval==D0centre) ) )
		ustruct_files = []
		for ll in range(-Nwindow,Nwindow+1):
			for dd in range(-Nwindow,Nwindow+1):
				Lcurrent = cellsize[Lcentreidx + ll]
				D0current = dval[D0centreidx + dd]
				searchkey = '{}/mc_sides*_diam{}um_fpert0.1_npert*_dval{}um2ms.conf_0.traj'.format(trajdir,format(Lcurrent,'.1f'),format(D0current,'.2f'))
				ustruct_buffer = gb.glob(searchkey)
				ustruct_files = ustruct_files + ustruct_buffer	
		
		# Create element of processing list
		processinglist.append([bvalseq,deltaseq,Deltaseq,Nspins,Nsteps,Tdur,ustruct_files,ii])     # Append each slice list and create a longer list of MRI slices whose processing will run in parallel
		# Print how long it took
		toc = time.time()
		print('                         ii = {} / {} in {} sec'.format(ii+1,len(ustruct_list),toc-tic))
	print('')


	print('++++++++ analysing {} microstructures for {}'.format(len(ustruct_list),outstring))
	print('')

	outlist = []
	for ii in range(0,len(ustruct_list)):
		# Measure how long it takes
		tic = time.time()
		# Compute MRI signals
		synresults = syn.getSignals(processinglist[ii])
		outlist.append(synresults)
		# Measure how long it took
		toc = time.time()
		print('                         ii = {} in {} sec'.format(ii,toc-tic))
	print('')


	# Save
	print('++++++++ saving outputs')
	print('')
	hfile = open('{}.sig.icell.SNRinf.bin'.format(outstring),'wb')
	pk.dump(outlist,hfile) 
	hfile.close()

	np.savetxt('{}.bval'.format(outstring), [bvalseq], fmt='%.1f', delimiter=' ')
	np.savetxt('{}.gDur'.format(outstring), [deltaseq], fmt='%.1f', delimiter=' ')
	np.savetxt('{}.gSep'.format(outstring), [Deltaseq], fmt='%.1f', delimiter=' ')

	sys.exit(0)
	

