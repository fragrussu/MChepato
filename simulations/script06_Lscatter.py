## Plot relationship between (D,K) and L at different SNR

import pickle as pk
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

bmax = 1000.0
snr_vals = [np.inf,20.0]   # SNR to plot; select from np.inf, 80.0, 40.0, 20.0
gDur_vals = [20.0,40.0,20.0,10.0,20.0]
gSep_vals = [25.0,50.0,50.0,50.0,75.0]
tDiff_vals = [18.3, 26.7, 43.3, 46.7, 68.3]
plt_title = ['(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)']

## Loop through different b-values
plot_count = 1
for ss in range(0,len(snr_vals)):
	
	try:
		SNR = '{}'.format(int(snr_vals[ss]))
		SNRstr = '= {}'.format(int(snr_vals[ss]))
	except:
		SNR = '{}'.format(snr_vals[ss])
		SNRstr = '→ ∞'
		
	
	## Loop through different MR protocols
	for nn in range(0,len(gDur_vals)):

		delta = gDur_vals[nn]
		Delta = gSep_vals[nn]
		Tdiff = tDiff_vals[nn] 
		

		# Load (D,K)
		h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.AdcAkcS0.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
		dataDtKtS0 = pk.load(h)
		h.close()

		# Load (L,D0)
		h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.ustruct.bin'.format(delta,Delta,bmax,delta,Delta,bmax),'rb')
		dataLD0 = pk.load(h)
		h.close()
		
		# Reformat D0, Dt, Kt as arrays
		Nustruct = len(dataLD0)
		L = np.zeros(Nustruct)
		for qq in range(0,Nustruct):
			L[qq] = dataLD0[qq][0]
		Dt = dataDtKtS0[0]
		Kt = dataDtKtS0[1]

		
		# Scatter plot
		plt.subplot(2,5,plot_count)
		plt.scatter(Dt,Kt,s=2.5,c=L,vmin=14, vmax=56)
		plt.xlim(0.0,2.40)
		plt.ylim(-5.0,10.0)
		plt.yticks(ticks=[-4,-2,0,2,4,6,8,10],fontsize=12)
		plt.xticks(ticks=[0.0,0.4,0.8,1.2,1.6,2.0,2.4],fontsize=12)
		plt.xlabel('$D\,\,\,\,\,[\mu m^2/ms] $',labelpad=0,fontsize=15)
		plt.ylabel('$K$',labelpad=-5,fontsize=15)
		plt.title('{}'.format(plt_title[plot_count-1]),fontsize=18,fontweight='bold')
		cbar = plt.colorbar(pad=0.015,shrink=0.85)
		cbar.set_label('$[\mu m]$', rotation=90, labelpad=-35,fontsize=11)
		cbar.set_ticks([14,56])
		cbar.ax.set_title('$L$',fontsize=15)
		cbar.ax.tick_params(labelsize=13)
		plt.text(0.95, -3, 'Δ-δ/3 = {} ms\nδ = {} ms\nΔ = {} ms\nSNR {}'.format(Tdiff,delta,Delta,SNRstr), horizontalalignment='left',verticalalignment='center',fontsize=11,fontweight='bold',bbox=dict(boxstyle = 'square',facecolor = [0.85,0.85,0.85]))

		plot_count = plot_count + 1

# Set up figure
plt.subplots_adjust(left=0.038, bottom=0.06, right=0.99, top=0.96, wspace=0.15, hspace=0.25)
fig = plt.gcf()
fig.set_size_inches([17.0,9.0])

# Show
plt.savefig('script06_Lscatter_bmax{}_900dpi.png'.format(bmax),dpi=900)
plt.show()
