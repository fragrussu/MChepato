### Print ranges for ADC and AKC for different MR protocols
import pickle as pk
import numpy as np


# SNR level:
snrval = 'inf'    # Select either '20' or 'inf'

# List of protocols
protlist = []
protlist.append('./results_gDur10.0_gSep50.0_bmin100.0_bmax1000.0_Nb7/syn_d10.0_D50.0_bmin100.0_bmax1000.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur10.0_gSep50.0_bmin100.0_bmax1500.0_Nb7/syn_d10.0_D50.0_bmin100.0_bmax1500.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur10.0_gSep50.0_bmin100.0_bmax2000.0_Nb7/syn_d10.0_D50.0_bmin100.0_bmax2000.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur20.0_gSep25.0_bmin100.0_bmax1000.0_Nb7/syn_d20.0_D25.0_bmin100.0_bmax1000.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur20.0_gSep25.0_bmin100.0_bmax1500.0_Nb7/syn_d20.0_D25.0_bmin100.0_bmax1500.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur20.0_gSep25.0_bmin100.0_bmax2000.0_Nb7/syn_d20.0_D25.0_bmin100.0_bmax2000.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur20.0_gSep50.0_bmin100.0_bmax1000.0_Nb7/syn_d20.0_D50.0_bmin100.0_bmax1000.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur20.0_gSep50.0_bmin100.0_bmax1500.0_Nb7/syn_d20.0_D50.0_bmin100.0_bmax1500.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur20.0_gSep50.0_bmin100.0_bmax2000.0_Nb7/syn_d20.0_D50.0_bmin100.0_bmax2000.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur20.0_gSep75.0_bmin100.0_bmax1000.0_Nb7/syn_d20.0_D75.0_bmin100.0_bmax1000.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur20.0_gSep75.0_bmin100.0_bmax1500.0_Nb7/syn_d20.0_D75.0_bmin100.0_bmax1500.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur20.0_gSep75.0_bmin100.0_bmax2000.0_Nb7/syn_d20.0_D75.0_bmin100.0_bmax2000.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur40.0_gSep50.0_bmin100.0_bmax1000.0_Nb7/syn_d40.0_D50.0_bmin100.0_bmax1000.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur40.0_gSep50.0_bmin100.0_bmax1500.0_Nb7/syn_d40.0_D50.0_bmin100.0_bmax1500.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))
protlist.append('./results_gDur40.0_gSep50.0_bmin100.0_bmax2000.0_Nb7/syn_d40.0_D50.0_bmin100.0_bmax2000.0_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(snrval))

# Print infor
for nn in range(0,len(protlist)):

	# Load ADC and AKC for current protocol
	h = open(protlist[nn],'rb'); 
	data = pk.load(h); 
	h.close()

	# Get ADC and AKC
	adc = data[0]
	akc = data[1]
	
	# Print ranges
	print('+++ Protocol {}'.format(protlist[nn]))
	print('    Median and 95% range of ADC: {:.3f} [{:.3f}; {:.3f}] um2/ms'.format( np.median(adc), np.percentile(adc,2.5) , np.percentile(adc,100.0 - 2.5) ) )
	print('    Median and 95% range of AKC: {:.3f} [{:.3f}; {:.3f}]'.format( np.median(akc), np.percentile(akc,2.5) , np.percentile(akc,100.0 - 2.5) ) )
	print('')
	
print('')
