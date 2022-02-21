## Show examples of (D0,L) prediction with polynomial interpolation (PolyMap) on the test set

import pickle as pk
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

bmax = 1000.0	
delta = 20.0
Delta = 25.0
Tdiff = np.round(10*(Delta - delta/3.0))/10.0
SNR = 40
try:
	SNRstr = '= {}'.format(int(SNR))
	SNRstr_outfile = '{}'.format(int(SNR))
except:
	SNRstr = '→ ∞'
	SNRstr_outfile = 'inf'


plt_title = ['(A) Ground truth','(B) Ground truth','(C) Prediction','(D) Prediction']

# Load (D,K)
print('** Loading ./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.AdcAkcS0.test.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR))
print('')
h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.AdcAkcS0.test.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
dataDtKtS0 = pk.load(h)
h.close()
Dt = dataDtKtS0[0]
Kt = dataDtKtS0[1]		


		
# Load ground truth (L,D0)
print('** Loading ./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.LD0gt.test.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
print('')
h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.LD0gt.test.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
dataLD0gt = pk.load(h)
h.close()
Lgt = dataLD0gt[0]
D0gt = dataLD0gt[1]


# Load estimators
print('** Loading ./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.AdcAkc_to_L.train.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
print('')
h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.AdcAkc_to_L.train.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
Lestimator = pk.load(h)
h.close()


print('** Loading ./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.AdcAkc_to_L.train.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
print('')
h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.AdcAkc_to_D0.train.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
D0estimator = pk.load(h)
h.close()


# Load predicted (L,D0)
print('** Loading ./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.LD0pred.test.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
print('')
h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.LD0pred.test.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
dataLD0pred = pk.load(h)
h.close()
Lpred = dataLD0pred[0]
D0pred = dataLD0pred[1]


# Print estimators
print('')
print('++++++++++++++')
print('Coefficients of L(D,K) expansion ')
print('    ')
print(' L = a0 + a1 D + a2 K + a3 D^2  + a4 K^2 + a5 D K + a6 D^2 K + a7 D K^2 + a8 D^3 + a9 K^3')
print('')
print(Lestimator)
print('')

print('')
print('++++++++++++++')
print('Coefficients of D0(D,K) expansion ')
print('    ')
print(' D0 = a0 + a1 D + a2 K + a3 D^2  + a4 K^2 + a5 D K + a6 D^2 K + a7 D K^2 + a8 D^3 + a9 K^3')
print('')
print(D0estimator)
print('')


# Plot figures

plot_count = 1
plt.subplot(2,2,plot_count)
plt.scatter(Dt,Kt,s=3.5,c=D0gt,vmin=0.3, vmax=2.3)
plt.xlim(0.0,2.40)
plt.ylim(-5.0,10.0)
plt.yticks(ticks=[-4,-2,0,2,4,6,8,10],fontsize=12)
plt.xticks(ticks=[0.0,0.4,0.8,1.2,1.6,2.0,2.4],fontsize=12)
plt.xlabel('$D\,\,\,\,\,[\mu m^2/ms] $',labelpad=0,fontsize=14)
plt.ylabel('$K$',labelpad=-5,fontsize=14)
plt.title('{}'.format(plt_title[plot_count-1]),fontsize=15,fontweight='bold')
cbar = plt.colorbar(pad=0.025,shrink=0.85)
cbar.set_label('$[\mu m^2/ms]$', rotation=90, labelpad=-25,fontsize=11)
cbar.set_ticks([0.3,2.3])
cbar.ax.set_title('$D_0$',fontsize=14)
cbar.ax.tick_params(labelsize=13)
plt.text(1.05, -2.5, 'Δ-δ/3 = {} ms\nδ = {} ms\nΔ = {} ms\nSNR {}'.format(Tdiff,delta,Delta,SNRstr), horizontalalignment='left',verticalalignment='center',fontsize=11,fontweight='bold',bbox=dict(boxstyle = 'square',facecolor = [0.85,0.85,0.85]))



plot_count = 3
plt.subplot(2,2,plot_count)
plt.scatter(Dt,Kt,s=3.5,c=D0pred,vmin=0.3, vmax=2.3)
plt.xlim(0.0,2.40)
plt.ylim(-5.0,10.0)
plt.yticks(ticks=[-4,-2,0,2,4,6,8,10],fontsize=12)
plt.xticks(ticks=[0.0,0.4,0.8,1.2,1.6,2.0,2.4],fontsize=12)
plt.xlabel('$D\,\,\,\,\,[\mu m^2/ms] $',labelpad=0,fontsize=14)
plt.ylabel('$K$',labelpad=-5,fontsize=14)
plt.title('{}'.format(plt_title[plot_count-1]),fontsize=15,fontweight='bold')
cbar = plt.colorbar(pad=0.025,shrink=0.85)
cbar.set_label('$[\mu m^2/ms]$', rotation=90, labelpad=-25,fontsize=11)
cbar.set_ticks([0.3,2.3])
cbar.ax.set_title('$D_0$',fontsize=14)
cbar.ax.tick_params(labelsize=13)



plot_count = 2
plt.subplot(2,2,plot_count)
plt.scatter(Dt,Kt,s=3.5,c=Lgt,vmin=14, vmax=56)
plt.xlim(0.0,2.40)
plt.ylim(-5.0,10.0)
plt.yticks(ticks=[-4,-2,0,2,4,6,8,10],fontsize=12)
plt.xticks(ticks=[0.0,0.4,0.8,1.2,1.6,2.0,2.4],fontsize=12)
plt.xlabel('$D\,\,\,\,\,[\mu m^2/ms] $',labelpad=0,fontsize=14)
plt.ylabel('$K$',labelpad=-5,fontsize=14)
plt.title('{}'.format(plt_title[plot_count-1]),fontsize=15,fontweight='bold')
cbar = plt.colorbar(pad=0.025,shrink=0.85)
cbar.set_label('$[\mu m]$', rotation=90, labelpad=-20,fontsize=11)
cbar.set_ticks([14,56])
cbar.ax.set_title('$L$',fontsize=14)
cbar.ax.tick_params(labelsize=13)




plot_count = 4
plt.subplot(2,2,plot_count)
plt.scatter(Dt,Kt,s=3.5,c=Lpred,vmin=14, vmax=56)
plt.xlim(0.0,2.40)
plt.ylim(-5.0,10.0)
plt.yticks(ticks=[-4,-2,0,2,4,6,8,10],fontsize=12)
plt.xticks(ticks=[0.0,0.4,0.8,1.2,1.6,2.0,2.4],fontsize=12)
plt.xlabel('$D\,\,\,\,\,[\mu m^2/ms] $',labelpad=0,fontsize=14)
plt.ylabel('$K$',labelpad=-5,fontsize=14)
plt.title('{}'.format(plt_title[plot_count-1]),fontsize=15,fontweight='bold')
cbar = plt.colorbar(pad=0.025,shrink=0.85)
cbar.set_label('$[\mu m]$', rotation=90, labelpad=-20,fontsize=11)
cbar.set_ticks([14,56])
cbar.ax.set_title('$L$',fontsize=14)
cbar.ax.tick_params(labelsize=13)




# Set up figure
plt.subplots_adjust(left=0.08, bottom=0.065, right=1.0, top=0.96, wspace=0.15, hspace=0.25)
fig = plt.gcf()
fig.set_size_inches([7.7,7.7])

# Show
plt.savefig('script07_ShowScatterPrediction_PolyMap_d{}_D{}_bmax{}_SNR{}_900dpi.png'.format(delta,Delta,bmax,SNRstr_outfile),dpi=900)
plt.show()


