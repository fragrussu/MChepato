import argparse, os, sys
import numpy as np
import nibabel as nib
import pandas
import pickle as pk

def GetPatchWiseMaps(flist,H,W,Hp,Wp,Zp,out,SzMax=42.0):
	'''
	This program converts a CSV file storing information on cell detection from high-resolution optical imaging of
	stained histological sections to a patch-wise parametric maps for comparison with Magnetic Resonance Imaging (MRI). 
	Author: Francesco Grussu, Vall d Hebron Institute of Oncology (<fgrussu@vhio.net><francegrussu@gmail.com>). 
	Copyright (c) 2021 Vall d Hebron Institute of Oncology (VHIO), Barcelona, Spain. All rights reserved.
	
	INTERFACE
	GetPatchWiseMaps(csvlist,H,W,Hp,Wp,out,SzMax=44.0)
	
	- flist: list of paths to CSV files storing the cell segmentation information (multiple files may be required for very large histological images). 
	           This code expect CSV files containing a column with variable name "Cell_Area", 
	           where different cell areas from all detected cells are reported in um^2, as well as two columns with variable names
	           "Centroid_X_um" and "Centroid_Y_um" storing the position of cells (in um) along the X (horizontal, i.e. image width) and Y 
	           (vertical, i.e. image height) direction
	- H:       field-of-view along the vertical direction (i.e. image height, in um) of the source 
	           histological image on which cells where segmented
	- W:       field-of-view along the horizontal direction (i.e. image width, in um) of the source 
	           histological image on which cells where segmented
	- Hp:      heigth of the patches in um, along the vertical direction (i.e. along the image height), within which statistics 
	           of cell size will be calculated. It should match the resolution along the same spatial direction of the MRI scan 
	           to which histological information is to be compared
	- Wp:      width of the patches in um, along the horizontal direction (i.e. along the image width), within which statistics 
	           of cell size will be calculated. It should match the resolution along the same spatial direction of the MRI scan 
	           to which histological information is to be compared
	- Zp:      thickness of the MRI slice to which the 2D histology is to be compared to (used to create the NIFTI header)
	- out:     root name of output files. There will be 4 output NIFTI files, with the following string added to the root name: 
		   *_Lum.nii -> cell size index (CSI), in um, with CSI = (<L^7>/<L^3>)^1/4, where L is the size of cells within a patch; 
		   *_Cellsmm2.nii -> cellularity map in cells/mm2, defined as number_of_cells_within_patch/patch_area_in_mm2;
		   *_FCellPatch.nii -> intra-cellular patch fraction 
		   *_ODeosin.nii ->  mean optical density of eosin.
		   The script will also store 3 pickle binary python files: 
		   *_CellSizePatches.bin, storing a list G where G[i][j] lists the sizes of all cells within patch (i,j); 
		   files ending as *_Lum.npy, *_Cellsmm2.npy, *_FCellPatch.npy, and *_ODeosin.npy storing the same maps 
		   as in *_Lum.nii (cell size index (CSI) map, in um), *_Cellsmm2.nii (cellularity map in cells/mm2), 
		   *_FCellPatch.nii (intra-cellular patch fraction) and *_ODeosin.nii (mean optical density of eosin) 
		   as NumPy binaries 
	- SzMax:   maximum realistic cell size in um (default: 42 um; cells larger than this value will be ignored)   		  

	'''
	### Load CSV files and create an array of cell areas and cell positions
	for nn in range(0,len(flist)):
		mycsv = flist[nn]
		print('')
		print('     ... loading file {} and stacking cell sizes and cell positions ...'.format(mycsv))
		d = pandas.read_csv(mycsv)
		if nn==0:
			area_array = d.Cell_Area
			Wpos_array = d.Centroid_X_um
			Hpos_array = d.Centroid_Y_um
			OD_array = d.Cell_Eosin_OD_mean
		else:
			area_array = np.concatenate([area_array,d.Cell_Area])
			Wpos_array = np.concatenate([Wpos_array,d.Centroid_X_um]) 
			Hpos_array = np.concatenate([Hpos_array,d.Centroid_Y_um])
			OD_array = np.concatenate([OD_array,d.Cell_Eosin_OD_mean])
		del d
	# Normalise optical density so that it ranges between 0 and 1
	OD_array = ( OD_array - np.min(OD_array) )/( np.max(OD_array) - np.min(OD_array) )
	
	### Process all cells and fill out a patchwise map
	## Get number of cells
	Ncells = area_array.size
	
	## Create empty patch-wise maps
	NHp = np.ceil(H/Hp)    # Number of patches along the image height
	NWp = np.ceil(W/Wp)    # Number of patches along the image width
	Hcorr = NHp*Hp         # Corrected histological image height after ceiling
	Wcorr = NWp*Wp         # Corrected histological image width after ceiling
	print('')
	print('     ... creating patch information:')
	print('                                    - {} patches of size {} um along image width (size: {} um)'.format(NWp,Wp,Wcorr))
	print('                                    - {} patches of size {} um along image heigth (size: {} um)'.format(NHp,Hp,Hcorr))
	Lmap = np.zeros((int(NHp),int(NWp)))   # Allocate cell size map
	Cmap = np.zeros((int(NHp),int(NWp)))   # Allocate cellularity map
	Fmap = np.ones((int(NHp),int(NWp)))    # Allocate intra-cellular patch fraction map
	ODmap = np.ones((int(NHp),int(NWp)))  # Allocate optical density of eosin
	
	## Find center of patches
	Hp_edges = np.linspace(0,Hcorr,int(NHp)+1) 
	Hp_centres = Hp_edges[0:int(NHp)] + Hp/2.0   # Position of the centre of each patch
	Wp_edges = np.linspace(0,Wcorr,int(NWp)+1) 
	Wp_centres = Wp_edges[0:int(NWp)] + Wp/2.0   # Position of the centre of each patch
	
	## Assign a label to all patches
	labelmap = np.zeros((int(NHp),int(NWp)))     # Map storing the label of each patch
	hhmap = np.zeros((int(NHp),int(NWp)))        # Map storing the patch number along the image heigth
	wwmap = np.zeros((int(NHp),int(NWp)))        # Map storing the patch number along the image width
	patch_count = 1
	for hh in range(0,int(NHp)):
		for ww in range(0,int(NWp)):
			labelmap[hh,ww] = patch_count
			hhmap[hh,ww] = hh
			wwmap[hh,ww] = ww
			patch_count = patch_count + 1
	
	## Assign cells to patches
	print('')
	print('     ... assigning cells to patches ...')
	cell_label = np.zeros(Ncells)
	for cc in range(0,Ncells):
	
		# Get cell position
		cc_hpos = Hpos_array[cc]
		cc_wpos = Wpos_array[cc]
		
		# Find patch to which the cell belongs
		ww_found = np.argmin( np.abs(Wp_centres - cc_wpos) )
		hh_found = np.argmin( np.abs(Hp_centres - cc_hpos) )
		
		# Store the label of the patch 
		cell_label[cc] = labelmap[hh_found,ww_found]	
			
	### Process patches
	print('')
	print('     ... calculating patch-wise statistics ...')
	
	# Create list to store cell size in each patch
	cell_list = [[]]*int(NHp)
	for hh in range(0,int(NHp)):
		cell_list[hh] = [[]]*int(NWp)
	
	# Loop through patches
	for hh in range(0,int(NHp)):
		for ww in range(0,int(NWp)):
			
			# Get the label of current patch
			hh_ww_label = labelmap[hh,ww]
			
			# Get histology stats from current patch
			hh_ww_OD = OD_array[cell_label==hh_ww_label]                                  # Array of optical densities of eosin
			hh_ww_A = area_array[cell_label==hh_ww_label]                                 # Array of areas in um^2
			hh_ww_L = (2.0/np.sqrt(np.pi))*np.sqrt(hh_ww_A)                               # Array of sizes in um
			hh_ww_A[hh_ww_L>SzMax] = np.nan                                               # Remove unrealistically large cells in cell area array
			hh_ww_L[hh_ww_L>SzMax] = np.nan                                               # Remove unrealistically large cells in cell size array
			cell_list[hh][ww] = hh_ww_L                                                   # Store array of sizes in um for the current patch
			
			# Compute metrics: CSI in um
			hh_ww_Lmri = ( np.nanmean(hh_ww_L**7)/np.nanmean(hh_ww_L**3) )**(1/4)         # Get patch-wise CSI map in um
			
			# Compute metrics: cellularity in cells/mm2
			hh_ww_Ncells = np.sum(~np.isnan(hh_ww_A))                                     # Number of cells within patch
			if (np.isnan(hh_ww_Lmri)):
				hh_ww_cell = 0.0
			else:
				hh_ww_cell = hh_ww_Ncells/(1e-6*Hp*Wp)                                # Get cellularity in cells/mm2
			
			# Compute metrics: intra-cellular patch fraction
			hh_ww_iFrac = np.nansum(hh_ww_A)/(Hp*Wp)
			if (np.isnan(hh_ww_iFrac)):
				hh_ww_iFrac = 0.0
			if (hh_ww_iFrac>1.0):
				hh_ww_iFrac = 1.0
				
			# Compute metrics: optical mean optical density of eosin	
			hh_ww_ODmean = 1.0 - np.nansum( (1.0 - hh_ww_OD)*hh_ww_A)/np.nansum(hh_ww_A)  # Get mean optical density of eosin
			if(np.isnan(hh_ww_ODmean)):
				hh_ww_ODmean = 1.0
			
			# Save in the 2D patch-wise maps
			Lmap[hh,ww] = hh_ww_Lmri
			Cmap[hh,ww] = hh_ww_cell
			Fmap[hh,ww] = hh_ww_iFrac
			ODmap[hh,ww] = hh_ww_ODmean
	
	### Save output patch-wise maps
	print('')
	print('     ... saving output patch-wise maps ...')
	
	## As python binaries
	np.save('{}_Lum.npy'.format(out),Lmap)                    # Cell size index as Numpy binary
	np.save('{}_Cellsmm2.npy'.format(out),Cmap)               # Cellularity as Numpy binary
	np.save('{}_FCellPatch.npy'.format(out),Fmap)             # Intra-cellular fraction as Numpy binary
	np.save('{}_ODeosin.npy'.format(out),ODmap)               # Optical density of eosin
	h_file = open('{}_CellSizePatches.bin'.format(out),'wb')  # List of cell sizes for all patches 
	pk.dump(cell_list,h_file,pk.HIGHEST_PROTOCOL)             
	h_file.close()	
	
	## As NIFTIs
	# Create header affine
	Amat = np.eye(4)
	Amat[0,0] = Hp/1000.0
	Amat[1,1] = Wp/1000.0
	Amat[2,2] = Zp/1000.0

	# Create 3D data matrices and save them: cell size index
	Dmat = np.zeros((int(NHp),int(NWp),1))
	Dmat[:,:,0] = Lmap
	img = nib.Nifti1Image(Dmat, Amat)
	img.header.set_xyzt_units(2)
	img.set_data_dtype('float64')
	img.set_data_dtype('float64')
	nib.save(img,'{}_Lum.nii'.format(out))
		
	# Create 3D data matrices and save them: cellularity
	Dmat = np.zeros((int(NHp),int(NWp),1))
	Dmat[:,:,0] = Cmap
	img = nib.Nifti1Image(Dmat, Amat)
	img.header.set_xyzt_units(2)
	img.set_data_dtype('float64')
	img.set_data_dtype('float64')
	nib.save(img,'{}_Cellsmm2.nii'.format(out)) 
	
	# Create 3D data matrices and save them: intra-cellular patch fraction
	Dmat = np.zeros((int(NHp),int(NWp),1))
	Dmat[:,:,0] = Fmap
	img = nib.Nifti1Image(Dmat, Amat)
	img.header.set_xyzt_units(2)
	img.set_data_dtype('float64')
	img.set_data_dtype('float64')
	nib.save(img,'{}_FCellPatch.nii'.format(out))  
	
	# Create 3D data matrices and save them: eosin mean optical density
	Dmat = np.zeros((int(NHp),int(NWp),1))
	Dmat[:,:,0] = ODmap
	img = nib.Nifti1Image(Dmat, Amat)
	img.header.set_xyzt_units(2)
	img.set_data_dtype('float64')
	img.set_data_dtype('float64')
	nib.save(img,'{}_ODeosin.nii'.format(out))  

	

if __name__ == "__main__":

	
	### Print help and parse arguments
	parser = argparse.ArgumentParser(description='This program converts a CSV file storing information on cell detection from high-resolution optical imaging of stained histological sections to a patch-wise parametric maps for comparison with Magnetic Resonance Imaging (MRI). Author: Francesco Grussu, Vall d Hebron Institute of Oncology (<fgrussu@vhio.net><francegrussu@gmail.com>). Copyright (c) 2021 Vall d Hebron Institute of Oncology (VHIO), Barcelona, Spain. All rights reserved.')
	parser.add_argument('csv_list', help='paths to CSV files storing the cell segmentation information (multiple files may be required for very large histological images). In case of multiple files, these must be separated by a comma (,). This code expects a column with variable name "Cell_Area", where different cell areas from all detected cells are reported in um^2, as well as two columns with variable names "Centroid_X_um" and "Centroid_Y_um" storing the position of cells (in um) along the X (horizontal, i.e. image width) and Y (vertical, i.e. image height) direction')
	parser.add_argument('Hfov', help='field-of-view along the vertical direction (i.e. image height, in um) of the source histological image on which cells where segmented')
	parser.add_argument('Wfov', help='field-of-view along the horizontal direction (i.e. image width, in um) of the source histological image on which cells where segmented')
	parser.add_argument('Hpatch', help='height of the patches in um, along the vertical direction (i.e. along the image height), within which statistics of cell size will be calculated. It should match the resolution along the same spatial direction of the MRI scan to which histological information is to be compared')
	parser.add_argument('Wpatch', help='width of the patches in um, along the horizontal direction (i.e. along the image width), within which statistics of cell size will be calculated. It should match the resolution along the same spatial direction of the MRI scan to which histological information is to be compared')
	parser.add_argument('Zpatch', help='thickness of the MRI slice to which the 2D histology is to be compared to (used to create the NIFTI header)')
	parser.add_argument('out_base', help='root name of output files. There will be 4 output NIFTI files, with the following string added to the root name: *_Lum.nii -> cell size index (CSI), in um, with CSI = (<L^7>/<L^3>)^1/4, where L is the size of cells within a patch;  *_Cellsmm2.nii -> cellularity map in cells/mm2, defined as number_of_cells_within_patch/patch_area_in_mm2; *_FCellPatch.nii -> intra-cellular patch fraction *_ODeosin.nii ->  mean optical density of eosin. The script will also store 3 pickle binary python files: *_CellSizePatches.bin, storing a list G where G[i][j] lists the sizes of all cells within patch (i,j); files ending as *_Lum.npy, *_Cellsmm2.npy, *_FCellPatch.npy, and *_ODeosin.npy storing the same maps as in *_Lum.nii (cell size index (CSI) map, in um), *_Cellsmm2.nii (cellularity map in cells/mm2), *_FCellPatch.nii (intra-cellular patch fraction) and *_ODeosin.nii (mean optical density of eosin) as NumPy binaries')
	parser.add_argument('--szmax', metavar='<value>', default='42.0', help='maximum realistic cell size in um (default: 42 um; cells larger than this value will be ignored)')
	
	args = parser.parse_args()

	### Get input information
	instr = args.csv_list
	inlist = instr.split(',')
	Hfov = float(args.Hfov)
	Wfov = float(args.Wfov)
	Hpatch = float(args.Hpatch)
	Wpatch = float(args.Wpatch)
	Zpatch = float(args.Zpatch)
	szmax = float(args.szmax)
	out_base = args.out_base
	
	### Print feedback
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('                             getPatchMapFromQuPath.py                               ')
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	print('')
	print('')
	print('* List of CSV files to process: {}'.format(inlist))
	print('* Source histological image height: {} um'.format(Hfov))
	print('* Source histological image width: {} um'.format(Wfov))
	print('* Target patch height: {} um'.format(Hpatch))
	print('* Target patch width: {} um'.format(Wpatch))
	print('* Target MRI slice thickness: {} um'.format(Zpatch))
	print('* Maximum allowed cell size: {} um'.format(szmax))
	print('* Output NIFTI files: {}_Lum.nii, {}_Cellsmm3.nii, {}_FCellPatch.nii, {}_ODeosin.nii'.format(args.out_base, args.out_base, args.out_base, args.out_base))
	print('* Output binary files: {}_Lum.npy, {}_Cellsmm3.npy, {}_FCellPatch.npy, {}_ODeosin.npy, {}_CellSizePatches.bin'.format(args.out_base, args.out_base, args.out_base, args.out_base, args.out_base))
	print('')
	print('')
	
	### Run code
	GetPatchWiseMaps(inlist,Hfov,Wfov,Hpatch,Wpatch,Zpatch,out_base,SzMax=szmax)
		
	### Done
	print('')
	print('     ... Done')
	print('')



