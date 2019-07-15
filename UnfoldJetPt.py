import rootpy
#from root_numpy import hist2array
from rootpy.plotting import Canvas,Legend
from rootpy.plotting.utils import draw
from rootpy.io import root_open
import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt
import numpy as np
from matplotlib.colors import LogNorm
from numpy import linspace
import ROOT


ROOT.gROOT.LoadMacro("src/RooUnfoldResponse.cxx")


f = root_open('CF_pPb_MC_legotrain/legotrain_397_20170804-0029_LHCb4_fix_CF_pPb_MC_ptHardMerged.root', 'read')
measJets = f.Get('AliJJetJtTask/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(0))
trueJets = f.Get('AliJJetJtTask/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(2))
jetPtResponse = f.Get('AliJJetJtTask/AliJJetJtMCHistManager/JetPtCorr/JetPtCorrNFin{:02d}'.format(0))
measJets.MarkerColor = 2
trueJets.MarkerColor = 1

fig, (ax1, ax2,ax3) = plt.subplots(3,1,figsize = (5,8), sharex=True, sharey=False,gridspec_kw = {'height_ratios':[2, 1,2]})
for ax in (ax1,ax2,ax3):
  ax.set_xscale('log')
  ax.set_xlabel(r'$p_{T}\left[GeV\right]$') #Add x-axis labels for bottom row
ax1.set_yscale('log')
ax2.set_yscale('linear')
ax1.set_ylabel(r'$\frac{dN}{dp_{T}}$')
ax2.set_ylabel('Ratio Meas/True')
ax3.set_ylabel(r'$p_{T,obs}\left[GeV\right]$')


ax1.text(0.02,0.005,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n  Full jets\n' r'Anti-$k_T$, R=0.4' '\nJet Cone',fontsize=7) #Add text to second subfigure, first parameters are coordinates in the drawn scale/units
  
for jets,label in zip((measJets,trueJets),('Measured','True')):
  rplt.errorbar(jets,xerr=False,emptybins=False,axes=ax1,label=label)
ax1.legend(loc = 'best')

ax1.grid(True) #Draw grid 
ratio = measJets.Clone()
ratio.Divide(trueJets)
rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax2)
ax1.set_xlim([5,500]) #Set x-axis limits
ax1.set_ylim([5e-12,2e0]) #Set y-axis limits

rplt.hist2d(jetPtResponse,axes=ax3,label="Response",norm=LogNorm())
x1 =linspace(0,500,50)
y1 = linspace(0,500,50)
plt.plot(x1,y1,'-r')
ax3.set_xlim([5,500])
ax3.set_ylim([5,500])
ax3.set_yscale('log')

plt.tight_layout()
plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
#plt.savefig(name,format='pdf') #Save figure
plt.show() #Draw figure on screen
