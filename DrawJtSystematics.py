import logging
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
import rootpy
import defs
import re
import matplotlib
from matplotlib import container
from ROOT import TGraphErrors
from rootpy.plotting import Canvas
from rootpy.plotting import Legend
from rootpy.plotting import Graph
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import math
from dataset import *
import numpy as np

d = dict(
    jetalg = r'Anti-$k_\mathrm{T}$, R=0.4',
    jettype = 'Full Jets',
    system = r'ALICE p-Pb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' ,
    trigger = 'kINT7/kEMCEJE',
    cut = r'$\left| \eta_\mathrm{jet} \right| < 0.25$'
    
)

labelsize= 15
styles = [20,24,25,27,30]
colors = [1,2,3,4,7]
colorsBox = ['k','r','g','b']
Rebin = 2
boxwidth = 2.5

def main(): 
  
  JetPtBins = [5,10,20,30,40,60,80,100,150,500]
  jetPtBins2 = [(JetPtBins[i],JetPtBins[i+1]) for i in range(8)]

  Njets = 8
  

  
  
  topcomment = "_systematics_Triggered"
  
  finderName = ["Full jets, Anti-k_\mathrm{T} R = 0.4","Charged Jets, Anti-k_\mathrm{T} R = 0.4"]
  finderType = ["Full","Charged"]
  setTitle = ["0","1","2","3","4","5"]
  finderR = {4,4}
  
  iC = 0
  logx = 1
  doWeight = 1
  mSize = 0.3
  iS = 0

  f = root_open("errors_test.root", 'read')



  start = 4
  errGraph = [f.Get('JetConeJtWeightBinNFin{:02d}JetPt{:02d}_Systematics'.format(iS,i)) for i in range(start,Njets)]  #Get jT histograms from file an array
  hJtSignalGraph = [f.Get('JetConeJtWeightBinNFin{:02d}JetPt{:02d}_Statistics'.format(iS,i)) for i in range(start,Njets)]  #Get jT histograms from file an array
  print(hJtSignalGraph)
  
  ylow = 5e-6
  yhigh = 5e9
  xlow = 0.1
  xhigh = 10

      
  
  ax = plt.gca()
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])

  for ij,h,h_sys,pT in zip(range(start,Njets),hJtSignalGraph,errGraph,jetPtBins2[start:]):
    print(ij)
    print("JetPt {}".format(pT))
    power = 2*ij-8
    scale = pow(10,power)
    h=defs.grrScale(h, scale)
    h.SetMarkerColor(colors[ij-start])
    label = r'${:02d}\:\mathrm{{GeV}} < p_{{\mathrm{{T,jet}}}} < {:02d}\:\mathrm{{GeV}} \times 10^{}$'.format(pT[0],pT[1],power)

    rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,label=label,fmt='+') #Plot jT histogram, 
    errorboxes = []
    n = h_sys.GetN()
    xs = h_sys.GetX()
    ys = h_sys.GetY()
    xerrs = h_sys.GetEX()
    yerrs = h_sys.GetEY()
    for x in (xs,ys,xerrs,yerrs):
      x.SetSize(n)
    for x, y, xe, ye in zip(xs, ys, xerrs, yerrs):
      rect = Rectangle((x - xe, (y - ye)*scale), xe*2,ye*2*scale)
      errorboxes.append(rect)
    pc = PatchCollection(errorboxes, facecolor=colorsBox[ij-4], alpha=0.5,edgecolor='None')
    ax.add_collection(pc)
    #print("({} - {})/9.0 + {} is {}".format(yhigh,ylow,ylow,(yhigh-ylow)/9.0+ylow))
  
  ax.text(1.7,3e6, d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['cut'],
            fontsize = 10)
      
  handles, labels = ax.get_legend_handles_labels()
  handles = [container.ErrorbarContainer(h,has_xerr=False,has_yerr=True) if isinstance(h, container.ErrorbarContainer) else h for h in handles]
  ax.legend(handles,labels,loc = 'lower left',numpoints=1)
  #ax.legend(loc = 'lower left')
  ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
  ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
  ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
  ytitle = r'$\frac{1}{N_{jets}}\;\frac{1}{j_\mathrm{T}}\;\frac{\mathrm{d} N}{\mathrm{d} j_\mathrm{T}}$'
  xtitle = r'$j_\mathrm{T}\left(\mathrm{GeV}/c\right)$'
  ax.set_ylabel(ytitle,fontsize=16)
  ax.set_xlabel(xtitle,fontsize=16)

  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  ax.set_xscale('log')
  ax.set_yscale('log')
  plt.savefig("{}.pdf".format(file),format='pdf') #Save figure

  plt.show() #Draw figure on screen
  

   
if __name__ == "__main__": main()

    