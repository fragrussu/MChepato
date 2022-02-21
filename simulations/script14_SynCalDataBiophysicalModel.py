#### Simulation of MRI signals and ADC computation to estimate c0 and c1 for the biophysical model in Equation 5
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
import random
import pickle as pk
import AdcAkcFit as dfit

# Random walk information
trajdir='/nfs/Data/Disk2/clfgrussu/FGrussu/fgrussu/simulations/data/meshes/prism/perturbed_randomwalks'    # Directory where random walks are stored
Nspins = 1000 # number of spins
Ntime = 3000  # number of time steps
Nsteps = Ntime + 1  # number of simulation interations
Tdur = 0.140  # total simulation duration in sec
Lcell_vals = np.array([11.0,15.5,23.0,29.0])
D0cell_vals = np.array([1.30,1.60,1.90,2.20])
delta_vals = np.array([10.0, 15.0, 20.0, 25.0, 30.0])
Delta_vals = np.array([30.0, 40.0, 50.0, 60.0, 70.0])
bmax = 800.0
Nbval = 5
Nshapes = 15  # Number of cell shape to use at fixed size

# Generate list of cells and sequence parameters to process
Lcell_list,D0cell_list,delta_list,Delta_list = np.meshgrid(Lcell_vals,D0cell_vals,delta_vals,Delta_vals)
Lcell_list = Lcell_list.flatten()
D0cell_list = D0cell_list.flatten()
delta_list = delta_list.flatten()
Delta_list = Delta_list.flatten()
Nproc = Lcell_list.size

# Output arrays to be filled
s0_list = np.zeros(Nproc)   # Non-DW signal
adc_list = np.zeros(Nproc)   # Apparent diffusion coefficient


# Loop through different cell sizes and sequences
for ii in range(0,Nproc):


	# Get current cell and sequence parameters
	L_ii = Lcell_list[ii]
	D0_ii = D0cell_list[ii]
	gdur_ii = delta_list[ii]
	gsep_ii = Delta_list[ii]
	bvalseq = np.linspace(0.0,bmax,Nbval)
	deltaseq = 0.0*bvalseq + gdur_ii
	Deltaseq = 0.0*bvalseq + gsep_ii
	deltaseq[bvalseq==0.0] = 0.0
	Deltaseq[bvalseq==0.0] = 0.0

	# Find trajectories and use 10 out of 15 cells for signal synthesis
	searchkey = '{}/mc_sides*_diam{}um_fpert0.1_npert*_dval{}um2ms.conf_0.traj'.format(trajdir,format(L_ii,'.1f'),format(D0_ii,'.2f'))
	ustruct_files = gb.glob(searchkey)
	random.shuffle(ustruct_files)
	ustruct_files = ustruct_files[0:Nshapes]
	
	# Synthesise signals
	processing_info = ([bvalseq,deltaseq,Deltaseq,Nspins,Nsteps,Tdur,ustruct_files,0])
	tic1 = time.time()
	synresults = syn.getSignals(processing_info)    # Useful info: data[0] --> signal, data[1][0] --> b-vaues in s/mm2, data[1][1] --> grad dur delta in ms, data[1][2] --> grad sep Delta in ms 
	toc1 = time.time()
	
	# Fit ADC
	tic2 = time.time()
	dparams = dfit.adcfit(np.array(synresults[0]),bvalseq)
	toc2 = time.time()
	s0_list[ii] = dparams[0]
	adc_list[ii] = dparams[1]
	
	# Print info
	print('')
	print('** Iteration {}/{}:'.format(ii+1,Nproc))
	print('                * Cell/Sequence parameters')
	print('                cell size = {} um'.format(L_ii))
	print('                cell diffusivity = {} um2/ms'.format(D0_ii))
	print('                delta = {} ms'.format(gdur_ii))
	print('                Delta = {} ms'.format(gsep_ii))
	print('                * Results')
	print('                signal synthesis done in {} sec'.format(toc1-tic1))
	print('                fitting done in {} sec'.format(toc2-tic2))
	print('                s0 = {}'.format(dparams[0]))
	print('                adc = {} um2/ms'.format(dparams[1]))
	print('')
	print('')


# Save
output_storage = [Lcell_list,D0cell_list,delta_list,Delta_list,bvalseq,s0_list,adc_list]
hout = open('script14_NCellShapes{}_ResultList_L-D0-del-Del-bv-s0-adc.bin'.format(Nshapes),'wb')
pk.dump(output_storage,hout)
hout.close()



