## Show examples of (D0,L) prediction with both PolyMap and SigFit on the test set
import pickle as pk
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

bmax = 1000.0	
delta = 20.0
Delta = 25.0
Tdiff = np.round(10*(Delta - delta/3.0))/10.0
SNR = 20
try:
	SNRstr = '= {}'.format(int(SNR))
	SNRstr_outfile = '{}'.format(int(SNR))
except:
	SNRstr = '→ ∞'
	SNRstr_outfile = 'inf'


sig_type = 'ivim'   # Choose among 'icell' and 'ivim'
data_set = 'test'   # Choose among 'train' and 'test'

plt_title = ['(A) $D_0$ ground truth','(B) $D_0$ prediction ($PolyMap$)','(C) $D_0$ prediction ($SigFit$)','(D) $L$ ground truth','(E) $L$ prediction ($PolyMap$)','(F) $L$ prediction ($SigFit$)']

# Load (D,K)
h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.AdcAkcS0.{}.{}.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,data_set,sig_type,SNR),'rb')
dataDtKtS0 = pk.load(h)
h.close()
Dt = dataDtKtS0[0]
Kt = dataDtKtS0[1]		


		
# Load ground truth (L,D0)
h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.LD0gt.{}.{}.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,data_set,sig_type,SNR),'rb')
dataLD0gt = pk.load(h)
h.close()
Lgt = dataLD0gt[0]
D0gt = dataLD0gt[1]



# Load predicted (L,D0): polynomial mapping (PolyMap)
h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.LD0pred.{}.{}.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,data_set,sig_type,SNR),'rb')
dataLD0pred = pk.load(h)
h.close()
Lpred_poly = dataLD0pred[0]
D0pred_poly = dataLD0pred[1]


# Load predicted (L,D0): signal fitting (SigFit)
h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.SigFit.LD0pred.{}.{}.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,data_set,sig_type,SNR),'rb')
dataLD0pred = pk.load(h)
h.close()
Lpred_sigfit = dataLD0pred[0]
D0pred_sigfit = dataLD0pred[1]


# Plot figures
plot_count = 1
plt.subplot(2,3,plot_count)
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



plot_count = 2
plt.subplot(2,3,plot_count)
plt.scatter(Dt,Kt,s=3.5,c=D0pred_poly,vmin=0.3, vmax=2.3)
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



plot_count = 3
plt.subplot(2,3,plot_count)
plt.scatter(Dt,Kt,s=3.5,c=D0pred_sigfit,vmin=0.3, vmax=2.3)
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



plot_count = 4
plt.subplot(2,3,plot_count)
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
plt.text(1.05, -2.5, 'Δ-δ/3 = {} ms\nδ = {} ms\nΔ = {} ms\nSNR {}'.format(Tdiff,delta,Delta,SNRstr), horizontalalignment='left',verticalalignment='center',fontsize=11,fontweight='bold',bbox=dict(boxstyle = 'square',facecolor = [0.85,0.85,0.85]))


plot_count = 5
plt.subplot(2,3,plot_count)
plt.scatter(Dt,Kt,s=3.5,c=Lpred_poly,vmin=14, vmax=56)
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



plot_count = 6
plt.subplot(2,3,plot_count)
plt.scatter(Dt,Kt,s=3.5,c=Lpred_sigfit,vmin=14, vmax=56)
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
plt.subplots_adjust(left=0.06, bottom=0.065, right=0.98, top=0.96, wspace=0.15, hspace=0.25)
fig = plt.gcf()
fig.set_size_inches([11.55,7.7])

# Show
plt.savefig('script07_ShowScatterPrediction_PolyMap_SigFit_d{}_D{}_bmax{}_SNR{}_900dpi.png'.format(delta,Delta,bmax,SNRstr_outfile),dpi=900)
plt.show()

