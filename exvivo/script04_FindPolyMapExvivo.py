## Find a polynomial mapping (D,K) ---> (D0,L) for the specific MRI sequence used ex vivo on the 9.4T

# Load packages
import numpy as np
import pickle as pk
import AdcAkcFit as dfit
from matplotlib import pyplot as plt

## Input/Output information
signal_file = 'exvivomri_d10.0_D30.0/results_syn_d10.0_D30.0.sig.icell.SNRinf.bin'
ustr_file = 'exvivomri_d10.0_D30.0/results_syn_d10.0_D30.0.ustruct.bin'
bth = 1700  # Minimum b-value to be used for ADC and AKC calculation
kmin = -5.0     # Lower bound for kurtosis   
kmax = 10.0     # Upper bound for kurtosis
b0snr = 100     # SNR at b = 0; set it to np.inf if no noise should be added
myseed = 19102018   # Seed number for reproducibility
Lfit_file = 'exvivomri_d10.0_D30.0/LfitPoly_d10.0_D30.0.sig.icell.SNR{}.bmin{}.bin'.format(b0snr,bth)
D0fit_file = 'exvivomri_d10.0_D30.0/D0fitPoly_d10.0_D30.0.sig.icell.SNR{}.bmin{}.bin'.format(b0snr,bth)

## Load synthetic measurements and ground truth microstructural parameters
h = open(signal_file,'rb')
sall = pk.load(h)
h.close() 


h = open(ustr_file,'rb')
uall = pk.load(h)
h.close()


## Get b-values and number of microstructures
Nustr = len(sall)
bvals = sall[0][1][0]


## Fit ADC and AKC for all microstructures
D_array = np.zeros(Nustr)
K_array = np.zeros(Nustr)
D0_array = np.zeros(Nustr)
L_array = np.zeros(Nustr)

np.random.seed(myseed)
for uu in range(0,Nustr):

	# Get signals and b-values
	sig = sall[uu][0]
	
	# Add noise if SNR is not np.inf
	if(~np.isinf(b0snr)):
		sig = np.sqrt( (sig + (1.0/b0snr)*np.random.randn(sig.size))**2 + ((1.0/b0snr)*np.random.randn(sig.size))**2 )    # Add Rician noise if required
	
	# Select high b-values
	#sigth = np.copy(sig)
	#bvalsth = np.copy(bvals)
	sigth = sig[bvals>=bth]
	bvalsth = bvals[bvals>=bth]
	
	# Estimate ADC and Kurtosis
	outpars = dfit.akcfit(sigth, bvalsth)
	D_array[uu] = outpars[1]
	K_array[uu] = outpars[2]
	
	# Store microstructural parameters (D0,L) corresponding to ADC and Kurtosis
	L_array[uu] = uall[uu][0]
	D0_array[uu] = uall[uu][1]
	
	# Print information
	print('')
	print('** Microstructure {}/{}:'.format(uu+1,Nustr))
	print('                * Cell/Sequence parameters')
	print('                cell size = {} um'.format(uall[uu][0]))
	print('                cell diffusivity = {} um2/ms'.format(uall[uu][1]))
	print('                * Signal')
	print('                signal = {}'.format(sigth))
	print('                b-values = {}'.format(bvalsth))
	print('                b = 0 SNR = {}'.format(b0snr))
	print('                * Fitting')
	print('                s0 = {}'.format(outpars[0]))
	print('                adc = {} um2/ms'.format(outpars[1]))
	print('                akc = {}'.format(outpars[2]))
	print('')
	print('')
	


## Find polynomial mapping (D,K) ---> (D0,L)

# Remove observations where kurtosis hits the fitting boundaries
myidx_fil = np.ones(Nustr,dtype='int')
myidx_fil[K_array==kmin] = int(0)
myidx_fil[K_array==kmax] = int(0)
D_array_filt = D_array[myidx_fil==int(1)]
K_array_filt = K_array[myidx_fil==int(1)]
L_array_filt = L_array[myidx_fil==int(1)]
D0_array_filt = D0_array[myidx_fil==int(1)]
Ntrain_filt = K_array_filt.size
			
D_array_filt = np.reshape(D_array_filt,(Ntrain_filt,1))
K_array_filt = np.reshape(K_array_filt,(Ntrain_filt,1))
L_array_filt = np.reshape(L_array_filt,(Ntrain_filt,1))
D0_array_filt = np.reshape(D0_array_filt,(Ntrain_filt,1))
			
			
# Fit	
Qmat = np.concatenate(( 1.0 + 0*D_array_filt, D_array_filt, K_array_filt, D_array_filt**2, K_array_filt**2, D_array_filt*K_array_filt, (D_array_filt**2)*K_array_filt, (K_array_filt**2)*D_array_filt, D_array_filt**3,  K_array_filt**3 ),axis=1)
Lfit = np.matmul( np.linalg.pinv( Qmat ) , L_array_filt )
D0fit = np.matmul( np.linalg.pinv( Qmat ) , D0_array_filt )


# Predict
Lpred =  Lfit[0] + Lfit[1]*D_array_filt + Lfit[2]*K_array_filt + Lfit[3]*D_array_filt*D_array_filt + Lfit[4]*K_array_filt*K_array_filt + Lfit[5]*D_array_filt*K_array_filt + Lfit[6]*D_array_filt*D_array_filt*K_array_filt + Lfit[7]*D_array_filt*K_array_filt*K_array_filt + Lfit[8]*D_array_filt*D_array_filt*D_array_filt + Lfit[9]*K_array_filt*K_array_filt*K_array_filt

D0pred =  D0fit[0] + D0fit[1]*D_array_filt + D0fit[2]*K_array_filt + D0fit[3]*D_array_filt*D_array_filt + D0fit[4]*K_array_filt*K_array_filt + D0fit[5]*D_array_filt*K_array_filt + D0fit[6]*D_array_filt*D_array_filt*K_array_filt + D0fit[7]*D_array_filt*K_array_filt*K_array_filt + D0fit[8]*D_array_filt*D_array_filt*D_array_filt + D0fit[9]*K_array_filt*K_array_filt*K_array_filt




# Save results
print('')
print('                .... saving results')
print('')
	
h = open(Lfit_file,'wb')
pk.dump(Lfit,h)
h.close()

h = open(D0fit_file,'wb')
pk.dump(D0fit,h)
h.close()


# Plot figures
plot_count = 1
plt.subplot(2,2,plot_count)
plt.scatter(D_array_filt,K_array_filt,s=3.5,c=D0_array_filt,vmin=0.3, vmax=2.3)
plt.xlim(0.0,2.40)
plt.ylim(0.0,1.0)
plt.yticks(ticks=[0,0.25,0.5,0.75,1.0],fontsize=12)
plt.xticks(ticks=[0.0,0.4,0.8,1.2,1.6,2.0,2.4],fontsize=12)
plt.xlabel('$D\,\,\,\,\,[\mu m^2/ms] $',labelpad=0,fontsize=14)
plt.ylabel('$K$',labelpad=-5,fontsize=14)
plt.title('Ground truth',fontsize=15,fontweight='bold')
cbar = plt.colorbar(pad=0.025,shrink=0.85)
cbar.set_label('$[\mu m^2/ms]$', rotation=90, labelpad=-25,fontsize=11)
cbar.set_ticks([0.3,2.3])
cbar.ax.set_title('$D_0$',fontsize=14)
cbar.ax.tick_params(labelsize=13)




plot_count = 3
plt.subplot(2,2,plot_count)
plt.scatter(D_array_filt,K_array_filt,s=3.5,c=D0pred,vmin=0.3, vmax=2.3)
plt.xlim(0.0,2.40)
plt.ylim(0.0,1.0)
plt.yticks(ticks=[0,0.25,0.5,0.75,1.0],fontsize=12)
plt.xticks(ticks=[0.0,0.4,0.8,1.2,1.6,2.0,2.4],fontsize=12)
plt.xlabel('$D\,\,\,\,\,[\mu m^2/ms] $',labelpad=0,fontsize=14)
plt.ylabel('$K$',labelpad=-5,fontsize=14)
plt.title('Predicted',fontsize=15,fontweight='bold')
cbar = plt.colorbar(pad=0.025,shrink=0.85)
cbar.set_label('$[\mu m^2/ms]$', rotation=90, labelpad=-25,fontsize=11)
cbar.set_ticks([0.3,2.3])
cbar.ax.set_title('$D_0$',fontsize=14)
cbar.ax.tick_params(labelsize=13)



plot_count = 2
plt.subplot(2,2,plot_count)
plt.scatter(D_array_filt,K_array_filt,s=3.5,c=L_array_filt,vmin=14, vmax=56)
plt.xlim(0.0,2.40)
plt.ylim(0.0,1.0)
plt.yticks(ticks=[0,0.25,0.5,0.75,1.0],fontsize=12)
plt.xticks(ticks=[0.0,0.4,0.8,1.2,1.6,2.0,2.4],fontsize=12)
plt.xlabel('$D\,\,\,\,\,[\mu m^2/ms] $',labelpad=0,fontsize=14)
plt.ylabel('$K$',labelpad=-5,fontsize=14)
plt.title('Ground truth',fontsize=15,fontweight='bold')
cbar = plt.colorbar(pad=0.025,shrink=0.85)
cbar.set_label('$[\mu m]$', rotation=90, labelpad=-20,fontsize=11)
cbar.set_ticks([14,56])
cbar.ax.set_title('$L$',fontsize=14)
cbar.ax.tick_params(labelsize=13)




plot_count = 4
plt.subplot(2,2,plot_count)
plt.scatter(D_array_filt,K_array_filt,s=3.5,c=Lpred,vmin=14, vmax=56)
plt.xlim(0.0,2.40)
plt.ylim(0.0,1.0)
plt.yticks(ticks=[0,0.25,0.5,0.75,1.0],fontsize=12)
plt.xticks(ticks=[0.0,0.4,0.8,1.2,1.6,2.0,2.4],fontsize=12)
plt.xlabel('$D\,\,\,\,\,[\mu m^2/ms] $',labelpad=0,fontsize=14)
plt.ylabel('$K$',labelpad=-5,fontsize=14)
plt.title('Predicted',fontsize=15,fontweight='bold')
cbar = plt.colorbar(pad=0.025,shrink=0.85)
cbar.set_label('$[\mu m]$', rotation=90, labelpad=-20,fontsize=11)
cbar.set_ticks([14,56])
cbar.ax.set_title('$L$',fontsize=14)
cbar.ax.tick_params(labelsize=13)



# Set up figure
plt.subplots_adjust(left=0.08, bottom=0.065, right=1.0, top=0.96, wspace=0.15, hspace=0.25)
fig = plt.gcf()
fig.set_size_inches([7.7,7.7])
plt.show()


 
