### Estimation of c0 and c1 for the biophysical model in Equation 5
import numpy as np
import pickle as pk
from matplotlib import pyplot as plt
from sklearn import linear_model
import statsmodels.api as sm

## Load results
input_file = 'script14_NCellShapes1_ResultList_L-D0-del-Del-bv-s0-adc.bin'
print('')
print('*** Coefficient calculation with syntehtic data -- file {}'.format(input_file))
print('')
print('')
h = open(input_file,'rb')
data = pk.load(h)
h.close()
L = data[0]
D0 = data[1]
gdur = data[2]
gsep = data[3]
adc = data[6]
myvar0 = L*L*L*L/( D0*gdur*(gsep - gdur/3.0) ) 
myvar1 = -L*L*L*L*L*L/( D0*D0*gdur*gdur*(gsep - gdur/3.0) ) 



### Estimate coefficient via robust regression: model A, D = a0*L*L*L*L/( D0*gdur*(gsep - gdur/3.0) ) - a1*L*L*L*L*L*L/( D0*D0*gdur*gdur*(gsep - gdur/3.0) )
# Prepare data
xvar0_col = np.reshape(myvar0,(myvar0.size,1))
xvar1_col = np.reshape(myvar1,(myvar1.size,1))
xvar_mat = np.concatenate((xvar0_col,xvar1_col),axis=1)
yvar_col = np.reshape(adc,(adc.size,1))

# Regress
modelA = sm.OLS(yvar_col,xvar_mat).fit()
paramsA = modelA.params
confintA = modelA.conf_int()
mycoeff0_modelA = paramsA[0]
mycoeff1_modelA = paramsA[1]
mycoeff0qmin_modelA = confintA[0,0]
mycoeff0qmax_modelA = confintA[0,1]
mycoeff1qmin_modelA = confintA[1,0]
mycoeff1qmax_modelA = confintA[1,1]
adc_pred_modelA = mycoeff0_modelA*myvar0 + mycoeff1_modelA*myvar1

print('')
print('')
print('   -- Coefficient estimation:')
print('           c0 =  {} [{}; {}]'.format(mycoeff0_modelA,mycoeff0qmin_modelA,mycoeff0qmax_modelA))
print('           c1 =  {} [{}; {}]'.format(mycoeff1_modelA,mycoeff1qmin_modelA,mycoeff1qmax_modelA))
print('')
print('')




