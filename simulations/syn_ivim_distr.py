#### Library of functions used to synthesise diffusion MRI signals from random walks (averaging of 3 gradient directions)

import numpy as np
import multiprocessing
import sys
import os
import glob as gb
import pickle as pk
import matplotlib.pyplot as plt
import time


# Load random walks and store spin trajectories in a usable format
def loadTraj(Nspins,Nsteps,Tdur,filelist):

	# Create array of time
	T = np.linspace(0,Tdur,Nsteps) # array of time in sec

	# Allocate variable to store trajectories
	Nfiles = len(filelist)
	if(Nfiles==0):
		raise RuntimeError('I could not find any random walks! Error.')

	deltaXall = np.zeros((0,Nsteps))
	deltaYall = np.zeros((0,Nsteps))
	deltaZall = np.zeros((0,Nsteps))

	# Load all trajectories and stack files
	for ff in range(0,Nfiles):

		## load spin trajectories for current file
		data = np.array(np.fromfile(filelist[ff], dtype="float32"))    # Spin trajectories in mm
		data = 0.001*data # Spin trajectories in m

		## store trajectories in a matrix of size Nspins x Nsteps

		# x-component
		x = data[0:np.copy(data).size:3]
		X = np.zeros((Nspins,Nsteps))
		deltaX = np.zeros((Nspins,Nsteps))
		for tt in range(0,Nsteps):
		    X[:,tt] = x[tt:x.size:Nsteps]
		for nn in range(0,Nspins):
		    deltaX[nn,:] = X[nn,:] - X[nn,0]

		## y-component
		y = data[1:np.copy(data).size:3]
		Y = np.zeros((Nspins,Nsteps))
		deltaY = np.zeros((Nspins,Nsteps))
		for tt in range(0,Nsteps):
		    Y[:,tt] = y[tt:y.size:Nsteps]
		for nn in range(0,Nspins):
		    deltaY[nn,:] = Y[nn,:] - Y[nn,0]

		## z-component
		z = data[2:np.copy(data).size:3]
		Z = np.zeros((Nspins,Nsteps))
		deltaZ = np.zeros((Nspins,Nsteps))
		for tt in range(0,Nsteps):
		    Z[:,tt] = z[tt:z.size:Nsteps]
		for nn in range(0,Nspins):
		    deltaZ[nn,:] = Z[nn,:] - Z[nn,0]	

		## merge spins from current file with spins from previous files
		deltaXall = np.concatenate((deltaXall,deltaX),axis=0)
		deltaYall = np.concatenate((deltaYall,deltaY),axis=0)
		deltaZall = np.concatenate((deltaZall,deltaZ),axis=0)

	return deltaXall, deltaYall, deltaZall, T


# Simulate MRI signals from a set of random walks
def getSignals(inlist):


	# Get MRI sequence information
	bval = inlist[0]*1e6       # b-value: input in sec/mm2, converted to sec/m2    
	delta = inlist[1]*1e-3     # Gradient duration: input in msec, converted to sec
	Delta = inlist[2]*1e-3     # Gradient separation: input in msec, converted to sec

	# Get info for loading trajectories
	Nspinfile = inlist[3] # number of spins per files
	Nsteps = inlist[4]  # number of time steps
	Tdur = inlist[5]  # total simulation duration in sec
	filelist = inlist[6]  # List of files with MRI trajectories
	listid = inlist[7]   # Position in the batch input list


	# Calculate gradient strength
	gammar = 267.522187e6      # Gyromagnetic ratio in 1/(sec T)
	Nmeas = bval.size
	Gval = np.zeros(Nmeas)
	for mm in range(0,Nmeas):
		if(bval[mm]>0.0):
			Gval[mm] = np.sqrt( bval[mm] / ( gammar*gammar * delta[mm]*delta[mm] * (Delta[mm] - delta[mm]/3.0) ) )	

	# Simulate MRI signals
	mrisig = np.zeros(bval.shape)
	for mm in range(0,Nmeas):

		# b = 0 image
		if(bval[mm]==0.0):
			mrisig[mm] = 1.0
			
		# DW image
		else:
			
			# Synthesise signal for each cell in the input list individually
			nondwi_intracell = 0.0
			dwi_intracell = 0.0
			for qq in range(0,len(filelist)):	
			
				# Get exact cell size
				myfilepath = filelist[qq]
				myfilestr = os.path.split(myfilepath)[1]
				mycellsize = float(myfilestr[14:18])
				nondwi_intracell = nondwi_intracell + mycellsize*mycellsize*mycellsize   # Accumulate non-DWI signal from current cell to calculate signal from all cells
			
				# Load trajectories
				Xpos, Ypos, Zpos, T = loadTraj(Nspinfile,Nsteps,Tdur,[myfilepath])   # Random walk trajectories in meters for the current cell
				Nmolecules = Xpos.shape[0]   # Total number of water molecules
				dT = T[1] - T[0]   # duration of time step in sec			
				
				# Generate gradient waveform
				Gtime = np.zeros((1,T.size))
				Gtime[0,T<=delta[mm]] = Gval[mm]
				Gtime[0,T>Delta[mm]] = -Gval[mm]
				Gtime[0,T>Delta[mm]+delta[mm]] = -0.0

				# Allocate variables to store phase accruals for 3 orthogonal gradients
				phix = np.zeros(Nmolecules)
				phiy = np.zeros(Nmolecules)
				phiz = np.zeros(Nmolecules)

				# Calculate phase accruals and magnetic moments for gradient directions along x, y, z
				for nn in range(0,Nmolecules):
					phix[nn] = -1.0*gammar*dT*np.sum(Gtime[0,:]*Xpos[nn,:]) # phase acrual for nn-th spin when gradient direciton is [1 0 0]
					phiy[nn] = -1.0*gammar*dT*np.sum(Gtime[0,:]*Ypos[nn,:]) # phase acrual for nn-th spin when gradient direciton is [0 1 0]
					phiz[nn] = -1.0*gammar*dT*np.sum(Gtime[0,:]*Zpos[nn,:]) # phase acrual for nn-th spin when gradient direciton is [0 0 1]  				
		
				# Get complex magnetic moments for all spins
				momx = np.exp(-1j*phix)
				momy = np.exp(-1j*phiy)
				momz = np.exp(-1j*phiz)
				
				# Average over spin ensemble
				sx = np.mean(momx)
				sy = np.mean(momy)
				sz = np.mean(momz)
				
				# Simulate trace imaging
				strace = mycellsize*mycellsize*mycellsize*( np.abs(sx) + np.abs(sy) + np.abs(sz) )/3.0
				
				# Accumulate signal from current cell to calculate signal from all cells
				dwi_intracell = dwi_intracell + strace
			
			# Normalise intra-cellular signal by the non-DW signal (sum over cell volumes)
			dwi_intracell = dwi_intracell/nondwi_intracell
			
			# Store results
			bval_ms_per_um2 = bval[mm]*1e-9     # b-value in ms/um2, to be multiplied to diffuisivity in um2/ms
			mrisig[mm] = dwi_intracell          

		
		
	return [mrisig,inlist]   # Return a 2-element list with: 
	                         # mrisig --> MRI signals; 
	                         # inlist --> input information used to generate signals. Contains:
	                         #            inlist[0]: b-values in sec/mm2
	                         #            inlist[1]: gradient durations in msec
	                         #            inlist[2]: gradient separations in msec
	                         #            inlist[3]: number of spins per file
	                         #            inlist[4]: number of time steps
	                         #            inlist[5]: total simulation duration in msec
	                         #            inlist[6]: list of files with MRI trajectories
	                         #            inlist[7]: Position in the batch input list


