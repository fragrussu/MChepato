## Show prediction error on test set for polynomial interpolation (PolyMap)
import pickle as pk
import numpy as np
import matplotlib
from matplotlib import pyplot as plt

bmax = 1000.0
SNR = np.inf   # SNR to plot; select from np.inf, 80, 40, 20
gDur_vals = [20.0,40.0,20.0,10.0,20.0]
gSep_vals = [25.0,50.0,50.0,50.0,75.0]
tDiff_vals = [18.3, 26.7, 43.3, 46.7, 68.3]
plt_title = ['(A)','(B)','(C)','(D)','(E)','(F)','(G)','(H)','(I)','(J)']

## Loop through different b-values
plot_count = 1
try:
	SNRstr = '= {}'.format(int(SNR))
	SNRstrfile = '{}'.format(int(SNR))
except:
	SNRstr = '→ ∞'
	SNRstrfile = 'inf'
		
	
## Loop through different MR protocols
for nn in range(0,len(gDur_vals)):

	delta = gDur_vals[nn]
	Delta = gSep_vals[nn]
	Tdiff = tDiff_vals[nn] 
		

	# Load ground truth (L,D0) from validation/test set
	h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.LD0gt.test.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
	dataLD0_gt = pk.load(h)
	h.close()
	L_gt = dataLD0_gt[0]
	D0_gt = dataLD0_gt[1]
	
	# Load predicted (L,D0) from validation/test set
	h = open('./results_gDur{}_gSep{}_bmin100.0_bmax{}_Nb7/syn_d{}_D{}_bmin100.0_bmax{}_Nbval7.LD0pred.test.ivim.SNR{}.bin'.format(delta,Delta,bmax,delta,Delta,bmax,SNR),'rb')
	dataLD0_pred = pk.load(h)
	h.close()
	L = dataLD0_pred[0]
	D0 = dataLD0_pred[1]
	
	# Calculate errors
	Lerror = L - L_gt
	D0error = D0 - D0_gt
	
	# Calculate statistics for each discrete value of ground truth D0 and L
	D0_gt_vals = np.unique(D0_gt)
	L_gt_vals = np.unique(L_gt)
	
	D0_gt_vals_stat = np.zeros(len(D0_gt_vals))
	D0_gt_vals_up = np.zeros(len(D0_gt_vals))
	D0_gt_vals_down = np.zeros(len(D0_gt_vals))
	
	L_gt_vals_stat = np.zeros(len(L_gt_vals))
	L_gt_vals_up = np.zeros(len(L_gt_vals))
	L_gt_vals_down = np.zeros(len(L_gt_vals))
	
	for dd in range(0, len(D0_gt_vals)):
		vals = D0error[D0_gt==D0_gt_vals[dd]]	
		D0_gt_vals_stat[dd] = np.percentile(vals,50.0)
		D0_gt_vals_up[dd] = np.percentile(vals,75.0)
		D0_gt_vals_down[dd] = np.percentile(vals,25.0)	
	D0pol_x = np.concatenate((D0_gt_vals,np.flip(D0_gt_vals)))
	D0pol_y = np.concatenate((D0_gt_vals_up,np.flip(D0_gt_vals_down)))

	for dd in range(0, len(L_gt_vals)):
		vals = Lerror[L_gt==L_gt_vals[dd]]	
		L_gt_vals_stat[dd] = np.percentile(vals,50.0)
		L_gt_vals_up[dd] = np.percentile(vals,75.0)
		L_gt_vals_down[dd] = np.percentile(vals,25.0)	
	Lpol_x = np.concatenate((L_gt_vals,np.flip(L_gt_vals)))
	Lpol_y = np.concatenate((L_gt_vals_up,np.flip(L_gt_vals_down)))
	
	# Scatter plot: D0 error vs ground truth D0
	plt.subplot(2,5,nn+1)
	plt.scatter(D0_gt,D0error,s=10,c='black',marker='o',facecolor='none',label='All points for varying $L$')
	plt.plot(D0_gt_vals,D0_gt_vals_stat,'o',linewidth=4,markersize=6,color='orange',label='Median for varying $L$')
	plt.fill(D0pol_x,D0pol_y,alpha=0.25,facecolor='black',edgecolor='none',label='IQR for varying $L$')
	plt.xlim(0.2,2.4)
	plt.ylim(-1.2,1.2)
	plt.xticks(ticks=[0.3, 0.8, 1.3, 1.8, 2.3])
	plt.yticks(ticks=[-1.2 , -0.8, -0.4,  0.0,  0.4,  0.8, 1.0, 1.2  ])
	plt.xlabel('Ground truth  $D_0\,\,\,\,\,[\mu m^2/ms] $',labelpad=0,fontsize=13)
	plt.ylabel('$D_0$ error $[\mu m^2/ms]$',labelpad=0,fontsize=13)
	plt.title('{}'.format(plt_title[nn]),fontsize=16,fontweight='bold')
	plt.grid()
	if(nn==0):
		plt.legend()
	plt.text(0.3, -0.92, 'Δ-δ/3 = {} ms\nδ = {} ms\nΔ = {} ms\nSNR {}'.format(Tdiff,delta,Delta,SNRstr), horizontalalignment='left',verticalalignment='center',fontsize=11,fontweight='bold',bbox=dict(boxstyle = 'square',facecolor = [153.0/255.0,204.0/255.0,255.0/255.0]))

	
	# Scatter plot: L error vs ground truth L
	plt.subplot(2,5,nn+1+5)
	plt.scatter(L_gt,Lerror,s=10,c='black',marker='o',facecolor='none',label='All points for varying $D_0$')
	plt.plot(L_gt_vals,L_gt_vals_stat,'o',linewidth=4,markersize=6,color='orange',label='Median for varying $D_0$')
	plt.fill(Lpol_x,Lpol_y,alpha=0.25,facecolor='black',edgecolor='none',label='IQR for varying $D_0$')
	plt.xlim(13,57)
	plt.ylim(-30.0,20.0)
	plt.xticks(ticks=[14., 21., 28., 35., 42., 49., 56.])
	plt.yticks(ticks=[-30,-20,-10,-0,10,20])
	plt.xlabel('Ground truth  $L\,\,\,\,\,[\mu m] $',labelpad=0,fontsize=13)
	plt.ylabel('$L$ error $[\mu m]$',labelpad=0,fontsize=13)
	plt.title('{}'.format(plt_title[nn+5]),fontsize=16,fontweight='bold')
	plt.grid()
	if(nn==0):
		plt.legend()
		plt.text(15, -13.0, 'Δ-δ/3 = {} ms\nδ = {} ms\nΔ = {} ms\nSNR {}'.format(Tdiff,delta,Delta,SNRstr), horizontalalignment='left',verticalalignment='center',fontsize=11,fontweight='bold',bbox=dict(boxstyle = 'square',facecolor = [153.0/255.0,204.0/255.0,255.0/255.0]))
	else:
		plt.text(15, -23.0, 'Δ-δ/3 = {} ms\nδ = {} ms\nΔ = {} ms\nSNR {}'.format(Tdiff,delta,Delta,SNRstr), horizontalalignment='left',verticalalignment='center',fontsize=11,fontweight='bold',bbox=dict(boxstyle = 'square',facecolor = [153.0/255.0,204.0/255.0,255.0/255.0]))

	
	

# Set up figure
plt.subplots_adjust(left=0.05, bottom=0.06, right=0.99, top=0.96, wspace=0.30, hspace=0.30)
fig = plt.gcf()
fig.set_size_inches([17.0,9.0])

# Show
plt.savefig('script08_PredError_PolyMap.bmax{}.SNR{}.900dpi.png'.format(bmax,SNRstrfile),dpi=900)
plt.show()

