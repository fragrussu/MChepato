### Perform ordinal logistic regression for all protocols at fixed SNR to classify discretised cell sizes

import numpy as np
import pickle as pk
from statsmodels import api as sm
from sklearn.metrics import accuracy_score as getAcc
from sklearn.metrics import confusion_matrix as getCMat
from matplotlib import pyplot as plt

# Information on MRI protocols
gDur = [10.0,10.0,10.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,40.0,40.0,40.0]
gSep = [50.0,50.0,50.0,25.0,25.0,25.0,50.0,50.0,50.0,75.0,75.0,75.0,50.0,50.0,50.0]
bmaxvals = [1000.0,1500.0,2000.0,1000.0,1500.0,2000.0,1000.0,1500.0,2000.0,1000.0,1500.0,2000.0,1000.0,1500.0,2000.0]
bmin = 100.0
Nprot = len(gDur)
SNR = 20  # Choose among np.inf, 100, 80, 40, 20


# Information on cell sizes
Lmin = 14.0
Lmax = 56.0
Lbounds = np.linspace(Lmin,Lmax,4)    # We aim to classify cells as small, medium, large

# Kurtosis bounds
kmin = -5.0     # Lower bound for kurtosis   
kmax = 10.0     # Upper bound for kurtosis


# Number of permutation to obtain the null distribution of random classifications
Nperm = 1000

# Seed for reproducibility of random permutations
myseed = 19102018
np.random.seed(myseed)


# Loop through MR protocols
for nn in range(0,Nprot):

	# Get MR protocol info
	delta = gDur[nn]
	Delta = gSep[nn]
	bmax = bmaxvals[nn]


	### Training data
	# Load measured signal (D,K) of training set
	h = open('results_gDur{}_gSep{}_bmin{}_bmax{}_Nb7/syn_d{}_D{}_bmin{}_bmax{}_Nbval7.AdcAkcS0.train.ivim.SNR{}.bin'.format(delta,Delta,bmin,bmax,delta,Delta,bmin,bmax,SNR),'rb')
	dataDtKtS0_train = pk.load(h)
	h.close()
	D_train = dataDtKtS0_train[0]
	K_train = dataDtKtS0_train[1]
	Ntrain = D_train.size
	
	# Load ground truth cell size of training set
	h = open('results_gDur{}_gSep{}_bmin{}_bmax{}_Nb7/syn_d{}_D{}_bmin{}_bmax{}_Nbval7.LD0gt.train.ivim.SNR{}.bin'.format(delta,Delta,bmin,bmax,delta,Delta,bmin,bmax,SNR),'rb')
	dataLD0_train = pk.load(h)
	h.close()
	L_train = dataLD0_train[0]
	
	# Transform training cell size to a categorical variable
	Lcat_train = 2.0*np.ones(L_train.shape)
	Lcat_train[L_train<Lbounds[2]] = 1.0
	Lcat_train[L_train<Lbounds[1]] = 0.0
	
	# Training set: remove observations where kurtosis hits the fitting boundaries
	myidx_fil = np.ones(Ntrain,dtype='int')
	myidx_fil[K_train==kmin] = int(0)
	myidx_fil[K_train==kmax] = int(0)
	D_train_filt = D_train[myidx_fil==int(1)]
	K_train_filt = K_train[myidx_fil==int(1)]
	Lcat_train_filt = Lcat_train[myidx_fil==int(1)]
	Ntrain_filt = K_train_filt.size
	
	# Reshape all training data for fitting a multinomial logistic regression model
	D_train_filt = np.reshape(D_train_filt,(Ntrain_filt,1))
	K_train_filt = np.reshape(K_train_filt,(Ntrain_filt,1))
	Lcat_train_filt = np.reshape(Lcat_train_filt,(Ntrain_filt,1))
	Lcat_train_filt = np.squeeze(Lcat_train_filt)
	Qmat_train_filt = np.concatenate(( 1.0 + 0.0*D_train_filt, D_train_filt, K_train_filt, D_train_filt**2, K_train_filt**2, D_train_filt*K_train_filt, (D_train_filt**2)*K_train_filt, (K_train_filt**2)*D_train_filt, D_train_filt**3,  K_train_filt**3 ),axis=1)   # Design matrix for model fitting (training set)
	
	
	D_train = np.reshape(D_train,(Ntrain,1))
	K_train = np.reshape(K_train,(Ntrain,1))
	Lcat_train = np.reshape(Lcat_train,(Ntrain,1))
	Lcat_train = np.squeeze(Lcat_train)
	Qmat_train = np.concatenate(( 1.0 + 0.0*D_train, D_train, K_train, D_train**2, K_train**2, D_train*K_train, (D_train**2)*K_train, (K_train**2)*D_train, D_train**3,  K_train**3 ),axis=1)   # Design matrix for model fitting (training set)
	
	
	
	### Validation/test data
	# Load measured signal (D,K) of validation/test set
	h = open('results_gDur{}_gSep{}_bmin{}_bmax{}_Nb7/syn_d{}_D{}_bmin{}_bmax{}_Nbval7.AdcAkcS0.test.ivim.SNR{}.bin'.format(delta,Delta,bmin,bmax,delta,Delta,bmin,bmax,SNR),'rb')
	dataDtKtS0_test = pk.load(h)
	h.close()
	D_test = dataDtKtS0_test[0]
	K_test = dataDtKtS0_test[1]
	Ntest = D_test.size
	
	# Load ground truth cell size of validation/test set
	h = open('results_gDur{}_gSep{}_bmin{}_bmax{}_Nb7/syn_d{}_D{}_bmin{}_bmax{}_Nbval7.LD0gt.test.ivim.SNR{}.bin'.format(delta,Delta,bmin,bmax,delta,Delta,bmin,bmax,SNR),'rb')
	dataLD0_test = pk.load(h)
	h.close()
	L_test = dataLD0_test[0]
	
	# Transform validation/test cell size to a categorical variable
	Lcat_test = 2.0*np.ones(L_test.shape)
	Lcat_test[L_test<Lbounds[2]] = 1.0
	Lcat_test[L_test<Lbounds[1]] = 0.0
	
	
	# Reshape all validation/test data for prediction with a trained multinomial logistic regression model
	D_test = np.reshape(D_test,(Ntest,1))
	K_test = np.reshape(K_test,(Ntest,1))
	Lcat_test = np.reshape(Lcat_test,(Ntest,1))
	Lcat_test = np.squeeze(Lcat_test)
	Qmat_test = np.concatenate(( 1.0 + 0.0*D_test, D_test, K_test, D_test**2, K_test**2, D_test*K_test, (D_test**2)*K_test, (K_test**2)*D_test, D_test**3,  K_test**3 ),axis=1)   # Design matrix for prediction (test set)


	#### Fitting of a multinomial logistic regression model and prediction with the trained model	
	model = sm.MNLogit(Lcat_train_filt, Qmat_train_filt, missing='drop')
	res = model.fit(disp=False)   # Fit multinomial logistic model with verbose "off"
	pout_train = res.predict(Qmat_train)     # Probabilistic output: training set
	pout_test = res.predict(Qmat_test)       # Probabilistic output: test set
	Lcat_train_pred = np.zeros(Ntrain)
	Lcat_test_pred = np.zeros(Ntest)
	for qq in range(0,Ntrain):
		pout_val = pout_train[qq]
		Lcat_train_pred[qq] = float(np.argmax(pout_val))
	for qq in range(0,Ntest):
		pout_val = pout_test[qq]
		Lcat_test_pred[qq] = float(np.argmax(pout_val))
	del model
	del res
	
	### Computer classification figures for both training / validation set
	acc_train = getAcc(Lcat_train,Lcat_train_pred)
	cmat_train = getCMat(Lcat_train,Lcat_train_pred,normalize='all')
	acc_test = getAcc(Lcat_test,Lcat_test_pred)
	cmat_test = getCMat(Lcat_test,Lcat_test_pred,normalize='all')
	
	### Get confidence interval for accuracy on validation set due to chance via permutations
	acc_test_Rndvals = np.zeros(Nperm)
	for pp in range(0,Nperm):
		modelRnd = sm.MNLogit(np.random.permutation(Lcat_train), Qmat_train, missing='drop')
		resRnd = modelRnd.fit(disp=False)            # Fit multinomial logistic model with verbose "off"
		pout_testRnd = resRnd.predict(Qmat_test)     # Probabilistic output: test set
		Lcat_test_predRnd = np.zeros(Ntest)
		for qq in range(0,Ntest):
			pout_valRnd = pout_testRnd[qq]
			Lcat_test_predRnd[qq] = float(np.argmax(pout_valRnd))
		del modelRnd
		del resRnd
		acc_test_Rndvals[pp] = getAcc(Lcat_test,Lcat_test_predRnd)
	
	### Print classification results
	print('')
	print('++++ MRI protocol results_gDur{}_gSep{}_bmin{}_bmax{}_Nb7/syn_d{}_D{}_bmin{}_bmax{}_Nbval7.LD0gt.test.ivim.SNR{}.bin'.format(delta,Delta,bmin,bmax,delta,Delta,bmin,bmax,SNR)) 
	print('     gDur = {} ms, gSep = {} ms, bmin = {} s/mm2, bmax = {} s/mm2, SNR = {}'.format(delta,Delta,bmin,bmax,SNR))
	print('         * Accuracy on training set: {}'.format(acc_train))
	print('         * Accuracy on test set: {}'.format(acc_test))
	print('         * 95% CI random class.: [{}; {}]'.format( np.percentile(acc_test_Rndvals,2.5) , np.percentile(acc_test_Rndvals,100.0 - 2.5) ))
	print('         * Accuracy on test set class 0 (small  cells): {}'.format(cmat_test[0,0]/np.sum(cmat_test[0,:])))
	print('         * Accuracy on test set class 1 (medium cells): {}'.format(cmat_test[1,1]/np.sum(cmat_test[1,:])))
	print('         * Accuracy on test set class 2 (large  cells): {}'.format(cmat_test[2,2]/np.sum(cmat_test[2,:])))
	print('         * Full test set confusion matrix:')
	print(cmat_test)
	print('')
	print('')
	print('')


