import numpy as np
from scipy.optimize import minimize


def AdcAkcFobj(tissue,meas,bvals):
	sig = tissue[0]*np.exp(-bvals*tissue[1] + (1.0/6.0)*bvals*bvals*tissue[1]*tissue[1]*tissue[2])
	fobj = np.sum( (meas - sig)*(meas - sig) )
	return fobj
	
def AdcFobj(tissue,meas,bvals):
	sig = tissue[0]*np.exp(-bvals*tissue[1] )
	fobj = np.sum( (meas - sig)*(meas - sig) )
	return fobj


def voxelfit(mymeas,bval):
	'''
		pout = voxelfit(mymeas,bval)
		
		Fits S(b) = S0 exp( -b D + (1/6)*b^2*D^2*K ) to multi b-value measurements
		
		INPUT ARGUMENTS
		* mymeas: array storing signal measurements at different b-values
		* bval: array storing the corresponding b-values in s/mm2
		
		OUTPUT
		* pout: 3-element array, s.t. pout[0] = S0, pout[1] = D, pout[2] = K
			Range for parameters are 
			0  <= S0 <= 1
			0  <= D <= 2.4
			-5 <= K <= 10
		
		Author: Francesco Grussu <fgrussu@vhio.net>
			06 July 2021
	
	'''
	# Reformat MRI sequence information 
	Nmeas = bval.size          # Number of measurements
	Bcol = np.reshape(bval/1000.0,(Nmeas,1))  # b-values as column array expressed in ms/um2
	Bmat = np.concatenate((1.0 + 0.0*Bcol,-1.0*Bcol,Bcol*Bcol/6.0),axis=1)   # Sequence parameter matrix for linear parameter estimation
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
	pout = np.array([S0val,Dtval,Ktval])
	
	# Return output
	return pout



def akcfit(mymeas,bval):
	'''
		pout = akcfit(mymeas,bval)
		
		Fits S(b) = S0 exp( -b D + (1/6)*b^2*D^2*K ) to multi b-value measurements
		
		INPUT ARGUMENTS
		* mymeas: array storing signal measurements at different b-values
		* bval: array storing the corresponding b-values in s/mm2
		
		OUTPUT
		* pout: 3-element array, s.t. pout[0] = S0, pout[1] = D, pout[2] = K
			Range for parameters are 
			0  <= S0 <= 1
			0  <= D <= 2.4
			-5 <= K <= 10
		
		Author: Francesco Grussu <fgrussu@vhio.net>
			06 July 2021
	
	'''


	myout = voxelfit(mymeas,bval)
	return myout 




def adcfit(mymeas,bval):
	'''
		pout = adcfit(mymeas,bval)
		
		Fits S(b) = S0 exp( -b D ) to multi b-value measurements
		
		INPUT ARGUMENTS
		* mymeas: array storing signal measurements at different b-values
		* bval: array storing the corresponding b-values in s/mm2
		
		OUTPUT
		* pout: 2-element array, s.t. pout[0] = S0, pout[1] = D
			Range for parameters are 
			0  <= S0 <= 1
			0  <= D <= 2.4
		
		Author: Francesco Grussu <fgrussu@vhio.net>
			06 July 2021
	
	'''
	# Reformat MRI sequence information 
	Nmeas = bval.size          # Number of measurements
	Bcol = np.reshape(bval/1000.0,(Nmeas,1))  # b-values as column array expressed in ms/um2
	Bmat = np.concatenate((1.0 + 0.0*Bcol,-1.0*Bcol),axis=1)   # Sequence parameter matrix for linear parameter estimation
	Acol = np.reshape(np.log(mymeas),(Nmeas,1))   # Log-signals arranged as a column for matrix multiplication

	# Perform linear fitting
	params = np.matmul( np.linalg.pinv( Bmat ) , Acol )    # Linear parameter estimation
	S0guess = np.exp(params[0])          # Apparent proton density S0
	Dguess = params[1]                 # Apparent diffusion coefficient

	# Remove NaNs or Infs
	if(np.isinf(Dguess)):
		if(np.isneginf(Dguess)):
			Dguess = 0.01
		else:
			Dguess = 2.40
	if(np.isinf(S0guess)):
		if(np.isneginf(S0guess)):
			S0guess = 0.0
		else:
			S0guess = 1.0
	if(np.isnan(Dguess)):
		Dguess = 1.205   # Mid point of the fitting range
	if(np.isnan(S0guess)):
		S0guess = 0.5	# Mid point of the fitting range
		
	if(Dguess<0.0):
		Dguess = 0.0
	if(Dguess>2.4):
		Dguess = 2.4
		
	if(S0guess<0.0):
		S0guess = 0.0
	if(S0guess>1.0):
		S0guess = 1.0	
		
	# Perform non-linear fitting using linear fit as starting point			
	pguess = [S0guess,Dguess]
	pbounds=[(0.0,1.0),(0.0,2.4)]
	modelfit = minimize(AdcFobj,pguess,args=tuple([np.squeeze(mymeas),np.squeeze(bval/1000.0)]),bounds=pbounds)
	fit_exit = modelfit.success
	if fit_exit==True:
		pfit = modelfit.x
		S0val = pfit[0]
		Dtval = pfit[1]
	else:
		S0val = S0guess
		Dtval = Dguess
		
	# Store results to output array
	pout = np.array([S0val,Dtval])
	
	# Return output
	return pout



