#### Add IVIM contamination to noise free signals

import pickle as pk
import numpy as np
import glob as gb
import os

# Seed for reproducibility
myseed = 19102018
np.random.seed(myseed)

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


# Parameters for IVIM contamination
Divim_min = 15.0    # um2/ms
Divim_max = 60.0    # um2/ms

fivim_min = 0.05
fivim_max = 0.50

# Generate distribution of IVIM parameters
fv_array = fivim_min + (fivim_max - fivim_min)*np.random.rand(nustr)
Dv_array = Divim_min + (Divim_max - Divim_min)*np.random.rand(nustr)

# Loop through microstructures
print('')
for pp in range(0,len(dirlist)):

	# Print information
	print('   ++++  Sequence parameter set {}/{}'.format(pp+1,len(dirlist)))
	print('') 

	# Find signals
	myfile = gb.glob('{}/*.sig.icell.SNRinf.bin'.format(dirlist[pp]))
	myfile = os.path.splitext(myfile[0])
	myfile = os.path.splitext(myfile[0])
	myfile = os.path.splitext(myfile[0])
	myfile = os.path.splitext(myfile[0])
	myfile = os.path.basename(myfile[0])

	# Load signal structure
	print('++++++++   Processing MRI protocol {}/{}: {}'.format(pp+1,len(dirlist),dirlist[pp]))
	print('')
	print('      .... load signals ...')
	hfile = open('{}/{}.sig.icell.SNRinf.bin'.format(dirlist[pp],myfile),'rb')
	sig = pk.load(hfile)
	hfile.close()
	
	# Loop through microstructures
	ivim_list = []
	print('')
	print('      .... add contamination ...')
	for nn in range(0,nustr):
	
		# Get actual signals and b-values
		tissue_sig = sig[nn][0]
		bvals = sig[nn][1][0]   # this gives b-values in s/mm2
		bvals = bvals/1000.0    # this gives b-values in ms/um2
		
		# Synthesise IVIM contamination
		ivim_sig = np.exp(-1.0*bvals*Dv_array[nn])
		
		# Add contamination
		tot_sig = fv_array[nn]*ivim_sig + (1.0 - fv_array[nn])*tissue_sig
		
		# Store
		sig[nn][0] = tot_sig
		ivim_list.append([fv_array[nn],Dv_array[nn]])
		
		
	# Save contaminated signals and list of IVIM contamination parameters
	print('')
	print('      .... save ...')
	hfile = open('{}/{}.sig.ivim.SNRinf.bin'.format(dirlist[pp],myfile),'wb')
	pk.dump(sig,hfile) 
	hfile.close()
	
	hfile = open('{}/{}.vascpar.ivim.SNRinf.bin'.format(dirlist[pp],myfile),'wb')
	pk.dump(ivim_list,hfile) 
	hfile.close()
	
	# Done
	print('')
	
	
	
