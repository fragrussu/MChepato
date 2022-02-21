###### Add noise to noise-free signals (true intra-cellular and intra-cellular + IVIM contamination)

import pickle as pk
import numpy as np
import glob as gb
import os
from matplotlib import pyplot as plt

# Seed for reproducibility
myseed = 19102018
np.random.seed(myseed)

# SNR levels of interest
snrvals = np.array([np.inf, 100, 80, 40, 20])

# Number of total microstructures to study
Nustr = 1189

# Number of microstructures to be used for training
Ntrain = 700
Ntest = Nustr - Ntrain   # Number of microstructures to be used for testing
allidx = np.arange(Nustr)
np.random.shuffle(allidx)
trainidx = allidx[0:Ntrain]      # Microstructures to be used for finding the functions D0(D,K) and L(D,K)
testidx = allidx[Ntrain:Nustr]   # Microstructures where the functions D0(D,K) and L(D,K) will be tested
kmin = -5.0     # Lower bound for kurtosis   
kmax = 10.0     # Upper bound for kurtosis

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


# Loop through signal types: intra-cellular or IVIM-contaminated
for stype in ['icell','ivim']:


	# Print information
	print(' #####   Processing signal type {} out of icell (no. 1) and ivim (no. 2)'.format(stype))
	print('')
	print('') 
 
	# Loop through microstructures
	print('')
	for pp in range(0,len(dirlist)):

		# Print information
		print('   +++   Processing MRI protocol {}/{}: {}'.format(pp+1,len(dirlist),dirlist[pp]))
		print('')
		print('') 

		# Loop through SNR levels
		for ss in range(0,snrvals.size):
		
			print('')
			print('     *** SNR = {} ({}/{})'.format(snrvals[ss],ss+1,snrvals.size))
			print('')

			if(np.isinf(snrvals[ss])):
				snr_str = '{}'.format(snrvals[ss])
			else:
				snr_str = '{}'.format(int(snrvals[ss]))
			
			### Find data
			myfile = gb.glob('{}/*.sig.{}.SNR{}.bin'.format(dirlist[pp],stype,snr_str))
						
			myfile = os.path.splitext(myfile[0])
			myfile = os.path.splitext(myfile[0])
			myfile = os.path.splitext(myfile[0])
			myfile = os.path.splitext(myfile[0])
			myfile = os.path.basename(myfile[0])

			### Load ADC-AKC parameters
			print('')
			print('                .... load ADC-AKC ...')
			hfile = open('{}/{}.AdcAkcS0.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'rb')
			mripar = pk.load(hfile)
			hfile.close()
			myDt_all = mripar[0]      # ADC estimates for the entire data set (training + testing/validation)
			myKt_all = mripar[1]      # AKC estimates for the entire data set (training + testing/validation)
			myS0_all = mripar[2]      # S0 estimates for the entire data set (training + testing/validation)
			
			
			### Load microstructures and create arrays of cell size and cell diffusivity
			print('')
			print('                .... load microstructures ...')
			hfile = open('{}/{}.ustruct.bin'.format(dirlist[pp],myfile),'rb')
			ustr = pk.load(hfile)
			hfile.close()
			Lgt_all = np.zeros(Nustr)
			D0gt_all = np.zeros(Nustr)
			for qq in range(0,Nustr):
				Lgt_all[qq] = ustr[qq][0]
				D0gt_all[qq] = ustr[qq][1]
			
			### Estimate (Dt,Kt) --> (D0,L) functions mapping microstructure via linear fitting of polynomial functions
			print('')
			print('                .... estimating (D,K) --> (D0, L) mapping ...')
			# Extract training data points
			myDt_train = myDt_all[trainidx]
			myKt_train = myKt_all[trainidx]
			Lgt_train = Lgt_all[trainidx]
			D0gt_train = D0gt_all[trainidx]
			
			# Training set: remove observations where kurtosis hits the fitting boundaries
			myidx_fil = np.ones(Ntrain,dtype='int')
			myidx_fil[myKt_train==kmin] = int(0)
			myidx_fil[myKt_train==kmax] = int(0)
			myDt_train_filt = myDt_train[myidx_fil==int(1)]
			myKt_train_filt = myKt_train[myidx_fil==int(1)]
			Lgt_train_filt = Lgt_train[myidx_fil==int(1)]
			D0gt_train_filt = D0gt_train[myidx_fil==int(1)]
			Ntrain_filt = myKt_train_filt.size
			
			myDt_train_filt = np.reshape(myDt_train_filt,(Ntrain_filt,1))
			myKt_train_filt = np.reshape(myKt_train_filt,(Ntrain_filt,1))
			Lgt_train_filt = np.reshape(Lgt_train_filt,(Ntrain_filt,1))
			D0gt_train_filt = np.reshape(D0gt_train_filt,(Ntrain_filt,1))
			
			# Fit	
			Qmat = np.concatenate(( 1.0 + 0*myDt_train_filt, myDt_train_filt, myKt_train_filt, myDt_train_filt**2, myKt_train_filt**2, myDt_train_filt*myKt_train_filt, (myDt_train_filt**2)*myKt_train_filt, (myKt_train_filt**2)*myDt_train_filt, myDt_train_filt**3,  myKt_train_filt**3 ),axis=1)
			Lfit = np.matmul( np.linalg.pinv( Qmat ) , Lgt_train_filt )
			D0fit = np.matmul( np.linalg.pinv( Qmat ) , D0gt_train_filt )
			
			### Used the estimated (Dt,Kt) --> (D0,L) functions
			print('')
			print('                .... predicting using the (D,K) --> (D0, L) mapping')
			
			# Apply the estimated functions to the whole training set
			myLpred_train = Lfit[0] + Lfit[1]*myDt_all[trainidx] + Lfit[2]*myKt_all[trainidx] + Lfit[3]*myDt_all[trainidx]*myDt_all[trainidx] + Lfit[4]*myKt_all[trainidx]*myKt_all[trainidx] + Lfit[5]*myDt_all[trainidx]*myKt_all[trainidx] + Lfit[6]*myDt_all[trainidx]*myDt_all[trainidx]*myKt_all[trainidx] + Lfit[7]*myDt_all[trainidx]*myKt_all[trainidx]*myKt_all[trainidx] + Lfit[8]*myDt_all[trainidx]*myDt_all[trainidx]*myDt_all[trainidx] + Lfit[9]*myKt_all[trainidx]*myKt_all[trainidx]*myKt_all[trainidx]
			myD0pred_train = D0fit[0] + D0fit[1]*myDt_all[trainidx] + D0fit[2]*myKt_all[trainidx] + D0fit[3]*myDt_all[trainidx]*myDt_all[trainidx] + D0fit[4]*myKt_all[trainidx]*myKt_all[trainidx] + D0fit[5]*myDt_all[trainidx]*myKt_all[trainidx] + D0fit[6]*myDt_all[trainidx]*myDt_all[trainidx]*myKt_all[trainidx] + D0fit[7]*myDt_all[trainidx]*myKt_all[trainidx]*myKt_all[trainidx] + D0fit[8]*myDt_all[trainidx]*myDt_all[trainidx]*myDt_all[trainidx] + D0fit[9]*myKt_all[trainidx]*myKt_all[trainidx]*myKt_all[trainidx]
			
			# Apply the estimated functions to the whole test set
			myLpred_test = Lfit[0] + Lfit[1]*myDt_all[testidx] + Lfit[2]*myKt_all[testidx] + Lfit[3]*myDt_all[testidx]*myDt_all[testidx] + Lfit[4]*myKt_all[testidx]*myKt_all[testidx] + Lfit[5]*myDt_all[testidx]*myKt_all[testidx] + Lfit[6]*myDt_all[testidx]*myDt_all[testidx]*myKt_all[testidx] + Lfit[7]*myDt_all[testidx]*myKt_all[testidx]*myKt_all[testidx] + Lfit[8]*myDt_all[testidx]*myDt_all[testidx]*myDt_all[testidx] + Lfit[9]*myKt_all[testidx]*myKt_all[testidx]*myKt_all[testidx]
			myD0pred_test = D0fit[0] + D0fit[1]*myDt_all[testidx] + D0fit[2]*myKt_all[testidx] + D0fit[3]*myDt_all[testidx]*myDt_all[testidx] + D0fit[4]*myKt_all[testidx]*myKt_all[testidx] + D0fit[5]*myDt_all[testidx]*myKt_all[testidx] + D0fit[6]*myDt_all[testidx]*myDt_all[testidx]*myKt_all[testidx] + D0fit[7]*myDt_all[testidx]*myKt_all[testidx]*myKt_all[testidx] + D0fit[8]*myDt_all[testidx]*myDt_all[testidx]*myDt_all[testidx] + D0fit[9]*myKt_all[testidx]*myKt_all[testidx]*myKt_all[testidx]


			### Save results
			print('')
			print('                .... saving results')
			
			# Train
			h = open('{}/{}.AdcAkcS0.train.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'wb')
			pk.dump([myDt_all[trainidx],myKt_all[trainidx],myS0_all[trainidx]],h)
			h.close()
			
			h = open('{}/{}.AdcAkc_to_L.train.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'wb')
			pk.dump(Lfit,h)
			h.close()

			h = open('{}/{}.AdcAkc_to_D0.train.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'wb')
			pk.dump(D0fit,h)
			h.close()

			h = open('{}/{}.LD0gt.train.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'wb')
			pk.dump([Lgt_all[trainidx],D0gt_all[trainidx]],h)
			h.close()

			h = open('{}/{}.LD0pred.train.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'wb')
			pk.dump([myLpred_train,myD0pred_train],h)
			h.close()

			h = open('{}/{}.trainidx.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'wb')
			pk.dump(trainidx,h)
			h.close()
			
			# Test
			h = open('{}/{}.AdcAkcS0.test.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'wb')
			pk.dump([myDt_all[testidx],myKt_all[testidx],myS0_all[testidx]],h)
			h.close()

			h = open('{}/{}.LD0gt.test.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'wb')
			pk.dump([Lgt_all[testidx],D0gt_all[testidx]],h)
			h.close()

			h = open('{}/{}.LD0pred.test.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'wb')
			pk.dump([myLpred_test,myD0pred_test],h)
			h.close()

			h = open('{}/{}.testidx.{}.SNR{}.bin'.format(dirlist[pp],myfile,stype,snr_str),'wb')
			pk.dump(testidx,h)
			h.close()
			
			
		print('')
		print('')
		
		
	
	
