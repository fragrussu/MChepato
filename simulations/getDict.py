### Code used to synthesise signals according to the biophysical model of intra-cellular diffusion of Equation 5
import numpy as np


def syn(bval,gdur,gsep,nwords=50):

	# Constants to compute ADC
	a0 = 0.0013416743239695044
	a1 = 1.2593797077696043e-05
	
	# Convert b-values from s/mm2 to ms/um2
	bval = bval/1000.0
	
	## Create dictionary
	d0min = 0.3
	d0max = 2.3
	lmin = 14.0
	lmax = 56.0
	d0vals = np.linspace(d0min,d0max,nwords)
	lvals = np.linspace(lmin,lmax,nwords)	
	D0mat,Lmat = np.meshgrid(d0vals,lvals)
	D0array = D0mat.flatten()
	Larray = Lmat.flatten()
	Nmicro = D0array.size 

	## Generate signals
	Nmeas = bval.size
	sigdict = np.zeros((Nmicro,Nmeas))
	
	# Loop over microstructures
	for uu in range(0,Nmicro):
	
			
		# Get average D0 and L 
		lcell = Larray[uu]
		dcell = D0array[uu]
		
		for mm in range(0,Nmeas):
		
			delta = gdur[mm]
			Delta = gsep[mm]
			
			# Get ADC
			adcguess = a0*(lcell**4)/(dcell*delta*( Delta - delta/3.0 )) - a1*(lcell**6)/(dcell*dcell*delta*delta*( Delta - delta/3.0 ))   # Guess ADC with analytical formula
			if(adcguess>dcell):
				adcguess = dcell	     # Correct unplausible values that can be obtained at short diffusion times
			if(adcguess<0):
				adcguess = 0.0
			adc_mean = adcguess			
				
			
			# Get MRI signals for corresponding microstructure
			if(bval[mm]==0):
				sigdict[uu,mm] = 1.0
			else:
				sigdict[uu,mm] = np.exp(-bval[mm]*adc_mean)
				
	# Return synthetic signals and ground truth microstructures
	return sigdict, Larray, D0array
	
