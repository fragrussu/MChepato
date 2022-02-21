###### Get ADC and AKC for all synthetic microstructures, MRI protocols and SNR levels

import pickle as pk
import numpy as np
import glob as gb
import os
from scipy.optimize import minimize

### ADC-AKC decay objective function for non-linear fitting; tissue[0] = S0; tissue[1] = D; tissue[2] = K
def AdcAkcFobj(tissue,meas,bvals):
	sig = tissue[0]*np.exp(-bvals*tissue[1] + (1.0/6.0)*bvals*bvals*tissue[1]*tissue[1]*tissue[2])
	fobj = np.sum( (meas - sig)*(meas - sig) )
	return fobj

# Seed for reproducibility
myseed = 19102018
np.random.seed(myseed)

# SNR levels of interest
snrvals = np.array([np.inf, 100.0,80.0,40.0,20.0])

# Number of microstructures to study
nustr = 1189

# Folders to process
dirlist=[]
dirlist.append('results_gDur10.0_gSep50.0_bmin100.0_bmax1000.0_Nb7')
dirlist.append('results_gDur10.0_gSep50.0_bmin100.0_bmax1500.0_Nb7')
dirlist.append('results_gDur10.0_gSep50.0_bmin100.0_bmax2000.0_Nb7')
dirlist.append('results_gDur20.0_gSep25.0_bmin100.0_bmax1000.0_Nb7')
dirlist.append('results_gDur20.0_gSep25.0_bmin100.0_bmax1500.0_Nb7')
dirlist.append('results_gDur20.0_gSep25.0_bmin100.0_bmax2000.0_Nb7')
dirlist.append('results_gDur20.0_gSep50.0_bmin100.0_bmax1000.0_Nb7')
dirlist.append('results_gDur20.0_gSep50.0_bmin100.0_bmax1500.0_Nb7')
dirlist.append('results_gDur20.0_gSep50.0_bmin100.0_bmax2000.0_Nb7')
dirlist.append('results_gDur20.0_gSep75.0_bmin100.0_bmax1000.0_Nb7')
dirlist.append('results_gDur20.0_gSep75.0_bmin100.0_bmax1500.0_Nb7')
dirlist.append('results_gDur20.0_gSep75.0_bmin100.0_bmax2000.0_Nb7')
dirlist.append('results_gDur40.0_gSep50.0_bmin100.0_bmax1000.0_Nb7')
dirlist.append('results_gDur40.0_gSep50.0_bmin100.0_bmax1500.0_Nb7')
dirlist.append('results_gDur40.0_gSep50.0_bmin100.0_bmax2000.0_Nb7')



# Loop through microstructures
print('')
for pp in range(0,len(dirlist)):

	# Print information
	print('++++++++   Processing MRI protocol {}/{}: {}'.format(pp+1,len(dirlist),dirlist[pp]))
	print('')
	print('') 

	# Loop through SNR levels
	for ss in range(0,snrvals.size):
	
		# Print some info
		print('')
		print('  *** SNR = {} ({}/{})'.format(snrvals[ss],ss+1,snrvals.size))
		print('')

		# Loop: purely intra-cellular signals and signals with IVIM contamination
		for sigtype in ['icell','ivim']:
			
			# Print some info
			print('')
			print('      --- Type of signals: {}'.format(sigtype))
			print('')
			
			# Find signals
			try:
				myfile = gb.glob('{}/*.sig.{}.SNR{}.bin'.format(dirlist[pp],sigtype,int(snrvals[ss])))
			except:
				myfile = gb.glob('{}/*.sig.{}.SNR{}.bin'.format(dirlist[pp],sigtype,snrvals[ss]))			
			myfile = os.path.splitext(myfile[0])
			myfile = os.path.splitext(myfile[0])
			myfile = os.path.splitext(myfile[0])
			myfile = os.path.splitext(myfile[0])
			myfile = os.path.basename(myfile[0])

			# Load signal structure
			print('')
			print('              .... load signals ...')
			try:
				hfile = open('{}/{}.sig.{}.SNR{}.bin'.format(dirlist[pp],myfile,sigtype,int(snrvals[ss])),'rb')
			except:
				hfile = open('{}/{}.sig.{}.SNR{}.bin'.format(dirlist[pp],myfile,sigtype,snrvals[ss]),'rb')
			sig = pk.load(hfile)
			hfile.close()
			
			# Get b-values
			bval = sig[0][1][0]     # b-values in s/mm2
			Nmeas = bval.size          # Number of measurements
			Bcol = np.reshape(bval/1000.0,(Nmeas,1))  # b-values as column array expressed in ms/um2
			Bmat = np.concatenate((1.0 + 0.0*Bcol,-1.0*Bcol,Bcol*Bcol/6.0),axis=1)   # Sequence parameter matrix for linear parameter estimation
			
			# Loop through microstructures
			print('')
			print('              .... estimate ADC and AKC ...')
			S0val_array = np.zeros(nustr)
			Dtval_array = np.zeros(nustr)
			Ktval_array = np.zeros(nustr)
			mycounter = 1
			for uu in range(0,nustr):		
			
				mymeas = sig[uu][0]    # MRI signals
				Acol = np.reshape(np.log(mymeas),(Nmeas,1))   # Log-signals arranged as a column for matrix multiplication
				
				# Perform linear fitting
				params = np.matmul( np.linalg.pinv( Bmat ) , Acol )    # Linear parameter estimation
				S0guess = np.exp(params[0])          # Apparent proton density S0
				Dguess = params[1]                 # Apparent diffusion coefficient
				Kguess = params[2]/(Dguess*Dguess)   # Apparent kurtosis coefficient
				
				# Remove NaNs or Infs
				if(np.isinf(Dguess)):
					if(np.isneginf(Dguess)):
						Dguess = 0.01
					else:
						Dguess = 2.40
				if(np.isinf(Kguess)):
					if(np.isneginf(Kguess)):
						Kguess = -5.0
					else:
						Kguess = 10.0
				if(np.isinf(S0guess)):
					if(np.isneginf(S0guess)):
						S0guess = 0.0
					else:
						S0guess = 1.0
				if(np.isnan(Dguess)):
					Dguess = 1.205   # Mid point of the fitting range
				if(np.isnan(Kguess)):
					Kguess = 2.5	# Mid point of the fitting range
				if(np.isnan(S0guess)):
					S0guess = 0.5	# Mid point of the fitting range
					
				if(Dguess<0.0):
					Dguess = 0.0
				if(Dguess>2.4):
					Dguess = 2.4
					
				if(Kguess<-5.0):
					Kguess = -5.0
				if(Kguess>10.0):
					Kguess = 10.0
					
				if(S0guess<0.0):
					S0guess = 0.0
				if(S0guess>1.0):
					S0guess = 1.0	
					
				# Perform non-linear fitting using linear fit as starting point			
				pguess = [S0guess,Dguess,Kguess]
				pbounds=[(0.0,1.0),(0.0,2.4),(-5.0,10.0)]
				modelfit = minimize(AdcAkcFobj,pguess,args=tuple([np.squeeze(mymeas),np.squeeze(bval/1000.0)]),bounds=pbounds)
				fit_exit = modelfit.success
				if fit_exit==True:
					pfit = modelfit.x
					S0val = pfit[0]
					Dtval = pfit[1]
					Ktval = pfit[2]				
				else:
					S0val = S0guess
					Dtval = Dguess
					Ktval = Kguess
					
				# Store results to output array
				S0val_array[uu] = S0val
				Dtval_array[uu] = Dtval
				Ktval_array[uu] = Ktval
					
			# Save output list containing three arrays: outlist[0] --> ADC; outlist[1] --> AKC; outlist[2] --> S0
			print('')
			print('              .... save ...')
			print('')
			outlist = [Dtval_array,Ktval_array,S0val_array]
			try:
				hout = open('{}/{}.AdcAkcS0.{}.SNR{}.bin'.format(dirlist[pp],myfile,sigtype,int(snrvals[ss])),'wb')
			except:
				hout = open('{}/{}.AdcAkcS0.{}.SNR{}.bin'.format(dirlist[pp],myfile,sigtype,snrvals[ss]),'wb')
			pk.dump(outlist,hout) 
			hout.close()
			del S0val_array
			del Dtval_array
			del Ktval_array

		
				
	print('')
	print('')
	
	
	
	
