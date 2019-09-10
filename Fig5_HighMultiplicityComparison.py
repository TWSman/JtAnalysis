import logging
from matplotlib.pyplot import box
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
import rootpy
import defs
import re
import matplotlib
from matplotlib import container

from rootpy.plotting import Canvas
from rootpy.plotting import Legend
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import sys
from dataset import *
import ROOT
import math
labelsize = 15
mSize = 6
names = []
styles = [20,24,25,27,30]
colors = [1,2,3,4,7,6]
xlow = 0.1
colorsBox = ['k','r','g','b']

d = dict(
    jetalg = r'Anti-$k_\mathrm{T}$, R=0.4',
    jettype = 'Full Jets',
    system = r'ALICE p-Pb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' ,
    trigger = 'kINT7/kEMCEJE',
    cut = r'$\left| \eta_\mathrm{jet} \right| < 0.25$'
    
)
n_figs = 6
start= 4
end=8



def getBackgroundSystematic(signals,signalsRandom,divisor):
  """
  Creates a set of systematic error boxes for drawing
  
  Args:
    signals: Distributions from default background method
    signalsRandom: Distributions from alternative background method
    divisor: Distribution to be used for systematic error of ratio
  """
  systematics = []
  for h1,h2,h3 in zip(signals,signalsRandom,divisor):
    system = []
    N = h1.GetNbinsX()
    for i in range(1,N+1):
      x1 = h1.GetBinCenter(i)-h1.GetBinWidth(i)/2
      x2 = h1.GetBinCenter(i)+h1.GetBinWidth(i)/2
      y1 = h1.GetBinContent(i)
      y2 = h2.GetBinContent(i)
      y3 = h3.GetBinContent(i)
      
      error = math.fabs(y1-y2)
      ratio  =y1/y3
      ratioerror = error/y3
      system.append([x1,x2,y1-error,y1+error,y1,error,ratio,ratioerror])
    systematics.append(system)
  return systematics

def main(): 
  Rebin = 4
  #Load data
  Mixed_FullJets_R04_HM_01 = datasetMixed("V0A, 0 -   1%",NFIN=0,range=(1,5),filename="CF_pPb_legotrain/legotrain_CF_pPb_2274_20181219/legotrain_CF_pPb_2274_20181219_LHC13cde.root",directory='AliJJetJtTask_Central01/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE_Central01/AliJJetJtHistManager',color=colors[1],style=24,rebin=Rebin)
  Mixed_FullJets_R04_HM_10 = datasetMixed("V0A, 0 - 10%",NFIN=0, range=(1,5),filename="CF_pPb_legotrain/legotrain_CF_pPb_2305_20190109/legotrain_CF_pPb_2305_20190109_LHC13bcde.root",directory='AliJJetJtTask_Central10/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE_Central10/AliJJetJtHistManager',color=colors[2],style=24,rebin=Rebin)
  FullJets_R04_MB = dataset("Minimum Bias",NFIN=0, range=(1,8),filename="CF_pPb_legotrain/legotrain_CF_pPb_2749_20190822/legotrain_CF_pPb_2749_20190822_LHC13de.root",directory="AliJJetJtTask_kEMCEJE/AliJJetJtHistManager",color=colors[0],style=24,rebin=Rebin)
  FullJets_R04_HM_01_ZDC = dataset("ZDC, 0 -    1%",NFIN=0, range=(1,8),filename="CF_pPb_legotrain/legotrain_CF_pPb_2749_20190822/legotrain_CF_pPb_2749_20190822_LHC13de.root",directory="AliJJetJtTask_kEMCEJE_Central01/AliJJetJtHistManager",color=colors[4],style=24,rebin=Rebin)
  FullJets_R04_HM_10_ZDC = dataset("ZDC, 0 - 10%",NFIN=0, range=(1,8),filename="CF_pPb_legotrain/legotrain_CF_pPb_2749_20190822/legotrain_CF_pPb_2768_2019_0825LHC13d.root",directory="AliJJetJtTask_kEMCEJE_Central10/AliJJetJtHistManager",color=colors[5],style=24,rebin=Rebin)

  
  
  inclusive,jetPt = FullJets_R04_MB.getHist('JetConeJtWeightBin',jetpt = True)
  incs = [inclusive]
  datasets = [FullJets_R04_MB]
  datasets.append(Mixed_FullJets_R04_HM_10)
  datasets.append(FullJets_R04_HM_10_ZDC)


  for data in datasets[1:]:
    incs.append(data.getHist('JetConeJtWeightBin',jetpt = False))
  names = [data.name() for data in datasets]
  signals  = [data.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = False) for data in datasets]
  signals_randomBg  = [data.getSubtracted('JetConeJtWeightBin','BgRndmJtWeightBin',jetpt = False,randomBG=True) for data in datasets]

  systematics = [getBackgroundSystematic(h1, h2,signals[0]) for h1,h2 in zip(signals,signals_randomBg)]
  
  background = [data.getHist('BgJtWeightBin',jetpt=False,isBg=True) for data in datasets]



  drawSignal(signals,systematics,names,colors,styles,jetPt)

def drawSignal(signals,systematics,names,colors,styles,jetPt):

  if(n_figs == 2):
    fig,axs = defs.makeRatio(xlog=True,ylog=True,d=d,shareY=False,figsize = (5,7.5),grid=False)
  else:
    fig, axs = defs.makegrid(n_figs/2,2,xlog=True,ylog=True,d=d,shareY=False,figsize= (10,7.5) if n_figs == 4 else (n_figs*15/8,7.5) )
  axs = axs.reshape(n_figs)
  if(n_figs == 2):
    pT = jetPt[start]
    axs[0].text(0.8,7,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'] + '\n' + r'${:02d}\: < p_{{\mathrm{{T,jet}}}} < {:02d}\:\mathrm{{GeV}}/c$'.format(pT[0],pT[1]),fontsize = 11)
  else:
    axs[1].text(0.12,0.00002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'],fontsize = 11)
  for signal,system,name,color,j in zip(signals,systematics,names,colors,range(10)):
    for jT,syst,pT,ax,i in zip(signal[start:],system[start:],jetPt[start:],axs[0:n_figs/2],range(0,9)):
      jT.SetMarkerColor(color)
      jT.SetMarkerStyle(styles[j])
      jT.SetLineColor(color)
      plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=name,fmt='o',fillstyle='none') #Plot jT histogram, 
      line = plot.get_children()[0]
      line.set_markersize(mSize)
      if(styles[j] > 23):
        line.set_markerfacecolor('none')
        #line.set_markeredgecolor(color)
        line.set_color(color)
      if(n_figs > 2):
        ax.text(0.5,1e2,r'${:02d}\: < p_{{\mathrm{{T,jet}}}} < {:02d}\:\mathrm{{GeV}}/c$'.format(pT[0],pT[1])) 
      ax.set_xlim([xlow,10]) #Set x-axis limits
      ax.set_ylim([1e-5,2e3]) #Set y-axis limits
      errorboxes = []
      for box in syst:
        x1,x2,y1,y2,yc,error,ratio,ratioerror = box

        rect = Rectangle((x1,y1),x2-x1,y2-y1)
        errorboxes.append(rect)
      pc = PatchCollection(errorboxes, facecolor=colorsBox[j], alpha=0.5,edgecolor=colorsBox[j])
      ax.add_collection(pc)
      ax.grid(False)
          
    
    ratios = []
    if(j>0):
      for syst,sig,div,ax in zip(system[start:],signal[start:],signals[0],axs[n_figs/2:n_figs+1]):
        ratioBoxes = []
        N = div.GetNbinsX()
        for box,i in zip(syst,range(1,N)):
          y = div.GetBinContent(i)
          x1,x2,y1,y2,yc,error,ratio,ratioerror = box
          rect = Rectangle((x1,ratio-ratioerror),x2-x1,ratioerror*2)
          ratioBoxes.append(rect)
        
        pc = PatchCollection(ratioBoxes, facecolor=colorsBox[j], alpha=0.5,edgecolor=colorsBox[j])
        ax.add_collection(pc)
      
    if j > 0:
      for jT,div in zip(signal,signals[0]):
        h = jT.Clone()
        h.Divide(div)
        ratios.append(h)
      axs[n_figs/2].set_ylabel('Ratio to {}'.format(names[0]),fontsize=12) #Add y-axis labels to left- and righmost subfigures
      if(n_figs > 4): axs[-1].set_ylabel('Ratio to {}'.format(names[0]),fontsize=12)       
      for ratio,pT,ax in zip(ratios[start:],jetPt[start:],axs[n_figs/2:n_figs+1]):
        plot = rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram,
        ax.plot([0,10],[1,1],'k--')
     
        line = plot.get_children()[0]
        line.set_markersize(mSize)
        if(styles[j] > 23):
          line.set_markerfacecolor('none')
        line.set_color(color)  
        #if(i == 0):
        ax.set_yscale('linear')
        ax.set_xlim([xlow,10]) #Set x-axis limits
        ax.set_ylim([0.1,2.5]) #Set y-axis limits     
        ax.grid(False)
  handles, labels = axs[0].get_legend_handles_labels()
  handles = [container.ErrorbarContainer(h,has_xerr=False,has_yerr=True) if isinstance(h, container.ErrorbarContainer) else h for h in handles]
  axs[0].legend(handles,labels,loc = 'lower left',numpoints=1,prop={'family': 'monospace'})  
  #axs[0].legend(loc = 'lower left')
  axs[0].text(0.11,3e2,"ALICE",weight='bold')
  fig.align_labels()
  plt.savefig("PythonFigures/HighMJetConeJtSignalPtFrom{}To{}.pdf".format(start,end),format='pdf') #Save figure
  plt.show() #Draw figure on screen



  
if __name__ == "__main__": main()
