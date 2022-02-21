# Folder `exvivo`
The folder `exvivo` contains the following items:

* `script01_preprocess.sh`: script pre-processing the 9.4T MRI scans
* `script02_EstimateSNR.sh`: estimate SNR on the 9.4T MRI scans at b = 0, TE = 45 ms (results in `script02_EstimateSNR.RESULTS.txt`)
* `script03_synExVivoProtocol.sh`: script synthesising MRI measurements corresponding to the 9.4T MRI protocol, to estimate the PolyMap mapping (D,K) --> (D0,L); the  synthetic measurements are stored in ./exvivomri_d10.0_D30.0/results_syn_d10.0_D30.0.sig.icell.SNRinf.bin
* `script04_FindPolyMapExvivo.py`: script estimating the (D,K) ---> (D0,L) PolyMap polynomial mapping for the specific diffusion protocol and SNR used ex vivo
* `script05_PolyMap.py`: script computing (D0,L) in the fixed mouse livers with the PolyMap approach
* `script05_SigFit.py`: script computing (D0,L) in the fixed mouse livers with the SigFit approach
* `script06_segmentationQuPath.png`: screenshot of parameters used for automatic cell segmentation on histology through QuPath graphical interfaces
* `script07_TxtToCsv.sh`: script converting QuPath detection measurement files (.txt) to CSV (.csv)
* `script08_getPatchwiseHistomaps.sh`: script calculating patch-wise histological parametric maps for each of two MRI slices of two specimens
* `script09_warphisto2mri.sh`: script warping histological information to MRI
* `script10_MRIhisto_stats.sh`: script printing distributions of MRI and histological metrics in the two samples (results in `script10_MRIhisto_stats.RESULTS.txt`)
* `script11_ShowMaps.py`: script plotting with python Matplotlib histological and MRI parametric maps on top of a b = 0 MRI image for the two spacimens
* `script12_PlotFitQualityImages.py`: script plotting MRI images and PolyMap and SigFit image predictions given the fitted model parameters
* `script13_PlotFitQualitySignals.py`: script plotting MRI signal measurements for a representative voxel jointly with PolyMap and SigFit fitting
* `script14_SigFitFixD0.py`: SigFit estimation of cell size on the ex vivo MRI data performed when D0 is not estimated but fixed
* `script15_PlotCellSizeFixedD0.py`: plots cell size parametric maps obtained when D0 is not estimated but fixed
* `script16_SigFitFixD0_stats.sh`: script printing distributions of SigFit cell size estimates in the two samples obtained when D0 is not estimated but fixed (results in `script16_SigFitFixD0_stats.RESULTS.txt`)


The scripts listed above rely on the following additional tools and data:

* `synMRI_intracell_distrD0andL_snrinf_bvallist.py`: code synthesising MRI measurements from Monte Carlo random walks, used by `script03_synExVivoProtocol.sh`
* `exvivomri_d10.0_D30.0`: stores PolyMap (D,K) --> (D0,L) mappings obtained at SNR = 68 (median SNR for the WT case) and SNR = 100 (median SNR for the PDX case)
* `AdcAkcFit.py`: tools used for apparent diffusion/kurtosis coefficient fitting
* `getDict.py`: tools used for SigFit model fitting
* `getCSVfromQuPath.sh`: code converting a QuPath detection measurement file (.txt) to CSV (.csv) 
* `getPatchMapFromQuPath.py`: code computing patch-wise histological parametric maps from QuPath detection measurements in CSV format
* `warpPatchHisto2MRI.py`: code warping patch-wise histological parametric maps to MRI space using DiPy

