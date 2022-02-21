###### Add noise to noise-free signals (true intra-cellular and intra-cellular + IVIM contamination)

import pickle as pk
import numpy as np
import glob as gb
import os

# Seed for reproducibility
myseed = 19102018
np.random.seed(myseed)

# SNR levels of interest
snrvals = np.array([100.0,80.0,40.0,20.0])

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
	
		print('')
		print('  *** SNR = {} ({}/{})'.format(snrvals[ss],ss+1,snrvals.size))
		print('')

		### Signals with no IVIM contamination
		print('')
		print('      --- Intra-cellular signals')
		print('')
		
		# Find signals
		myfile = gb.glob('{}/*.sig.icell.SNRinf.bin'.format(dirlist[pp]))
		myfile = os.path.splitext(myfile[0])
		myfile = os.path.splitext(myfile[0])
		myfile = os.path.splitext(myfile[0])
		myfile = os.path.splitext(myfile[0])
		myfile = os.path.basename(myfile[0])

		# Load signal structure
		print('')
		print('              .... load signals ...')
		hfile = open('{}/{}.sig.icell.SNRinf.bin'.format(dirlist[pp],myfile),'rb')
		sig = pk.load(hfile)
		hfile.close()
		
		# Loop through microstructures
		print('')
		print('              .... add noise ...')
		for nn in range(0,nustr):
		
			# Get noise-free signals
			sig_nfree = sig[nn][0]
			
			# Add noise
			sig_noisy = np.sqrt( ( sig_nfree + (1.0/snrvals[ss])*np.random.randn(sig_nfree.size) )**2  +  ( (1.0/snrvals[ss])*np.random.randn(sig_nfree.size) )**2  )
			
			# Store
			sig[nn][0] = sig_noisy
			
			
		# Save noisy signals
		print('')
		print('              .... save ...')
		print('')
		hfile = open('{}/{}.sig.icell.SNR{}.bin'.format(dirlist[pp],myfile,int(snrvals[ss])),'wb')
		pk.dump(sig,hfile) 
		hfile.close()
			
		
		### Signals with IVIM contamination
		print('')
		print('      --- IVIM-contaminated signals')
		print('')
		
		# Find signals
		myfile = gb.glob('{}/*.sig.ivim.SNRinf.bin'.format(dirlist[pp]))
		myfile = os.path.splitext(myfile[0])
		myfile = os.path.splitext(myfile[0])
		myfile = os.path.splitext(myfile[0])
		myfile = os.path.splitext(myfile[0])
		myfile = os.path.basename(myfile[0])

		# Load signal structure
		print('')
		print('              .... load signals ...')
		hfile = open('{}/{}.sig.ivim.SNRinf.bin'.format(dirlist[pp],myfile),'rb')
		sig = pk.load(hfile)
		hfile.close()
		
		# Loop through microstructures
		print('')
		print('              .... add noise ...')
		for nn in range(0,nustr):
		
			# Get noise-free signals
			sig_nfree = sig[nn][0]
			
			# Add noise
			sig_noisy = np.sqrt( ( sig_nfree + (1.0/snrvals[ss])*np.random.randn(sig_nfree.size) )**2  +  ( (1.0/snrvals[ss])*np.random.randn(sig_nfree.size) )**2  )
			
			# Store
			sig[nn][0] = sig_noisy
			
			
		# Save noisy signals
		print('')
		print('              .... save ...')
		print('')
		hfile = open('{}/{}.sig.ivim.SNR{}.bin'.format(dirlist[pp],myfile,int(snrvals[ss])),'wb')
		pk.dump(sig,hfile) 
		hfile.close()
		
		
		# Done
		print('')
		print('')
		
	print('')
	print('')
	
	
	
	
