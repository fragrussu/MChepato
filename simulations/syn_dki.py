#### Library of functions used to synthesise diffusion MRI signals from random walks and for multiple gradient orientations (full DKI-like)

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
	bval = inlist[0]*1e6   # b-value: input in sec/mm2, converted to sec/m2    
	delta = inlist[1]*1e-3 # Gradient duration: input in msec, converted to sec
	Delta = inlist[2]*1e-3 # Gradient separation: input in msec, converted to sec
	grads = inlist[3]      # Matrix with gradient directions for each non-zero b-values (N gradients x 3, with grads[:,0] = x, grads[:,0] = y, grads[:,0] = z)
	Ngrads = grads.shape[0]   # Number of gradient directions
	
	# Get info for loading trajectories
	Nspinfile = inlist[4] # number of spins per files
	Nsteps = inlist[5]  # number of time steps
	Tdur = inlist[6]  # total simulation duration in sec
	filelist = inlist[7]  # List of files with MRI trajectories
	listid = inlist[8]   # Position in the batch input list

	# Calculate gradient strength
	gammar = 267.522187e6      # Gyromagnetic ratio in 1/(sec T)
	Nmeas = bval.size
	Gval = np.zeros(Nmeas)
	for mm in range(0,Nmeas):
		if(bval[mm]>0.0):
			Gval[mm] = np.sqrt( bval[mm] / ( gammar*gammar * delta[mm]*delta[mm] * (Delta[mm] - delta[mm]/3.0) ) )	

	# Create full scheme file-like MRI protocol information
	for mm in range(0,Nmeas):
		if(bval[mm]==0.0):
			bval_cur = np.array([0.0]) 
			Gval_cur = np.array([0.0])
			delta_cur = np.array([0.0])
			Delta_cur = np.array([0.0])
			gx_cur = np.array([0.0])
			gy_cur = np.array([0.0])
			gz_cur = np.array([0.0])
		else:
			bval_cur = bval[mm]*np.ones(Ngrads)
			Gval_cur = Gval[mm]*np.ones(Ngrads)
			delta_cur = delta[mm]*np.ones(Ngrads)
			Delta_cur = Delta[mm]*np.ones(Ngrads)
			gx_cur =  np.squeeze(grads[:,0])	
			gy_cur =  np.squeeze(grads[:,1])
			gz_cur =  np.squeeze(grads[:,2])
		
		try:
			bval_scheme = np.concatenate((bval_scheme,bval_cur))
			Gval_scheme = np.concatenate((Gval_scheme,Gval_cur))
			delta_scheme = np.concatenate((delta_scheme,delta_cur))
			Delta_scheme = np.concatenate((Delta_scheme,Delta_cur))
			gx_scheme = np.concatenate((gx_scheme,gx_cur)) 
			gy_scheme = np.concatenate((gy_scheme,gy_cur)) 
			gz_scheme = np.concatenate((gz_scheme,gz_cur)) 
			
		except:
			bval_scheme = np.copy(bval_cur)
			Gval_scheme = np.copy(Gval_cur)
			delta_scheme = np.copy(delta_cur)
			Delta_scheme = np.copy(Delta_cur)
			gx_scheme = np.copy(gx_cur)
			gy_scheme = np.copy(gy_cur)
			gz_scheme = np.copy(gz_cur)
			
	myscheme = [bval_scheme*1e-9,Gval_scheme,delta_scheme,Delta_scheme,gx_scheme,gy_scheme,gz_scheme]

	# Simulate MRI signals
	Ntot = bval_scheme.size
	mrisig = np.zeros(Ntot)
	for mm in range(0,Ntot):

		# b = 0 image
		if(bval_scheme[mm]==0.0):
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
				
				# Generate gradient waveform for all gradient components (x, y and z)
				Gtime = np.zeros((1,T.size))
				Gtime[0,T<=delta_scheme[mm]] = Gval_scheme[mm]
				Gtime[0,T>Delta_scheme[mm]] = -Gval_scheme[mm]
				Gtime[0,T>Delta_scheme[mm]+delta_scheme[mm]] = -0.0
				Gtimex = Gtime*gx_scheme[mm]
				Gtimey = Gtime*gy_scheme[mm]
				Gtimez = Gtime*gz_scheme[mm]
				
				# Variable to store phase accrual along gg-th direction
				phiarray = np.zeros(Nmolecules)
				
				# Calculate phase accruals and magnetic moments along gg-th gradient direction
				for nn in range(0,Nmolecules):
					phiarray[nn] = np.sum( -1.0*gammar*dT*( Gtimex[0,:]*Xpos[nn,:] + Gtimey[0,:]*Ypos[nn,:] + Gtimez[0,:]*Zpos[nn,:] ) ) 
				
				# Get complex magnetic moments for all spins
				magmom = np.exp(-1j*phiarray)
				
				# Average over spin ensemble
				scpx = np.mean(magmom)
				
				# Get volume-weighted MRI signal
				smag = mycellsize*mycellsize*mycellsize*( np.abs(scpx) )
				
				# Accumulate signal from current cell to calculate signal from all cells
				dwi_intracell = dwi_intracell + smag
			
			# Normalise intra-cellular signal by the non-DW signal (sum over cell volumes)
			dwi_intracell = dwi_intracell/nondwi_intracell
			
			# Store results
			mrisig[mm] = dwi_intracell          

		
	
		
	return [mrisig,myscheme]  # Return a 2-element list with: 
	                          # mrisig --> MRI signals; 
	                          # myscheme --> input information used to generate signals. Contains:
	                          #            myscheme[0]: b-values in ms/um2
	                          #            myscheme[1]: gradient strength in T/m
	                          #            myscheme[2]: gradient durations in sec
	                          #            myscheme[3]: gradient separations in sec     
	                          #            myscheme[4]: x components of the gradient
	                          #            myscheme[5]: y components of the gradient
	                          #            myscheme[6]: z components of the gradient



