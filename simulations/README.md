# Folder `simulations`
The folder `simulations` contains the following scripts used to perform the analyses presented in the paper:

* `runMC00_createvertices.py`: code used to generate the vertices of the meshes describing the synthetic cells, by perturbing regular triangularly-meshed prisms.
* `runMC01_createply.sh`: code used to create the .PLY files to be fed to the MCDC simulator
* `runMC02_runmc.sh`: code running the Monte Carlo diffusion simulation with MCDC (https://github.com/jonhrafe/MCDC_Simulator_public). It relies on three text files used to generate simulation configuration files (runMC02_runmc.confp1.conf, runMC2_runmc.confp2.conf, runMC2_runmc.confp3.conf)
* `script01run_gDur<...>_gSep<...>_bmin100.0_bmax<...>_Nb7.sh`: codes used to synthesise MRI signals for different clinically-feasible protocols
* `script02_addIVIM.py`: code adding intra-voxel incoherent motion (IVIM)-like contamination to intra-cellular signals
* `script03_addNoise.py`: code adding Rician noise to the synthetic MRI signal, after IVIM contamination
* `script04_getAdcAkc.py`: code estimating apparent diffusion/kurtosis coefficients (D,K) from MRI signals
* `script05_EstCellSize_PolyMap.py`: code establishing a polynomial mapping between (D,K) and cell diffusivity and cell size (D0,L) (PolyMap estimation)
* `script05_EstCellSize_SigFit.py`: code estimating cell diffusivity and cell size (D0,L) through direct fitting of a biophysical model of intra-cellular diffusion (SigFit estimation)
* `script06_D0andL_JointScatter.py`: code scattering (D,K) and colouring points according to D0 and L
* `script06_D0scatter.py`: code scattering (D,K) and colouring points according to D0
* `script06_Lscatter.py`: code scattering (D,K) and colouring points according to L
* `script07_ShowScatterPrediction_PolyMap.py`: code showing the PolyMap estimation (D,K) -> (D0,L) learnt on the training set at work on the validation set 
* `script07_ShowScatterPrediction_PolyMap_SigFit.py`: code showing PolyMap and SigFit estimation of (D0,L) by colouring (D,K) scatter plots according to ground truth and predicted (D0,L)
* `script08_PredError_PolyMap.py`: code scattering PolyMap (D0,L) prediction errors against ground truth (D0,L) on the validation set
* `script08_PredError_SigFit.py`: code scattering SigFit (D0,L) prediction errors against ground truth (D0,L) on the validation set (results stored in script09_OrdinalLogisticReg_CellSize.SNR<value>.RESULTS.txt)
* `script09_OrdinalLogisticReg_CellSize.py`: code performing multinomial logistic regression to classify a discrete number of cell sizes
* `script11a_run_<...>dir.sh`: code used to synthesise MRI signals for different gradient directions per b-value
* `script11b_analyse.py`: code used to fit ADC/AKC and full diffusion kurtosis tensor representations on one synthetic microstructural scenarios
* `script11c_plot.py`: code used to compare ADC/AKC and full diffusion kurtosis tensor fitting
* `script12a1_run.sh`: code used to synthesise MRI signals for 19 b-values for one microstructural scenario
* `script12a2_run.sh`: code used to synthesise MRI signals for 19 b-values for one microstructural scenario (different maximum b-value)
* `script12b_plot.py`: code used to visualise the impact of the number of b-values
* `script13_PrintAdcAkcRanges.py`: code used to extract ranges of apparent diffusion/kurtosis coefficients (results stored in script13_PrintAdcAkcRanges.RESULTS.SNR<value>.txt)
* `script14_SynCalDataBiophysicalModel.py`: code synthesising MRI signals and estimating intra-cellular ADC values to be used to estimate c0 and c1 from Equation 5 of the manuscript 
* `script15_FindConstantsBiophysicalModel.py`: estimation of c0 and c1 from Equation 5 of the manuscript 


The scripts listed above rely on the following additional tools:

* `AdcAkcFit.py`: tools used for apparent diffusion/kurtosis coefficient fitting
* `getDict.py`: tools used for SigFit model fitting
* `dir3.txt`: set of 3 isotropically-distributed gradient directions
* `dir9.txt`: set of 9 isotropically-distributed gradient directions, downloaded from http://www.emmanuelcaruyer.com/q-space-sampling.php
* `dir21.txt`: set of 21 isotropically-distributed gradient directions, downloaded from http://www.emmanuelcaruyer.com/q-space-sampling.php
* `dir30.txt`: set of 30 isotropically-distributed gradient directions, downloaded from http://www.emmanuelcaruyer.com/q-space-sampling.php
* `dir61.txt`: set of 61 isotropically-distributed gradient directions, downloaded from http://www.emmanuelcaruyer.com/q-space-sampling.php
* `runMC02_runmc.confp1.conf`: text file used for configuration on Monte Carlo simulations
* `runMC02_runmc.confp2.conf`: text file used for configuration on Monte Carlo simulations
* `runMC02_runmc.confp3.conf`: text file used for configuration on Monte Carlo simulations
* `script01_synMRI_intracell_distrD0andL_snrinf.py`: library of tools used to synthesise MRI signals for all microstructural scenarios by averaging 3 gradient directions
* `script11_synDKI_OneD0andLVal_SNRinf.py`: library of tools used to synthesise MRI signals for one microstructural scenario on several gradient directions, without averaging them
* `script12_synMRI_intracell_SubsetD0andL_snrinf.py`: library of tools used to synthesise MRI signals for one microstructural scenario by averaging 3 gradient directions
* `syn_dki.py`: tools used by `script11_synDKI_OneD0andLVal_SNRinf.py`
* `syn_ivim_distr.py`: tools used by `script12_synMRI_intracell_SubsetD0andL_snrinf.py`
* `script11b_analyse.adc.bin`: binary files used to perform plots via `script11c_plot.py`
* `script11b_analyse.akc.bin`: binary files used to perform plots via `script11c_plot.py`
* `script11b_analyse.md.bin`: binary files used to perform plots via `script11c_plot.py`
* `script11b_analyse.mk.bin`: binary files used to perform plots via `script11c_plot.py`
* `script14_NCellShapes1_ResultList_L-D0-del-Del-bv-s0-adc.bin`: binary files used to perform c0 and c1 estimation for Equation 5 with `script15_FindConstantsBiophysicalModel.py`

