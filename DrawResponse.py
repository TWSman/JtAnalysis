import rootpy
import defs
import re
import ROOT
import matplotlib
from rootpy.plotting import Canvas
from rootpy.plotting import Legend
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
from matplotlib.colors import PowerNorm
from matplotlib.colors import Normalize

import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
from dataset import *
from numpy import linspace
def main(): 
  print 'Number of arguments: ', len(sys.argv), 'arguments.'
  print 'Argument list:',str(sys.argv)
  filename = sys.argv[1]
  print "Input file: "
  print filename
  MC_FullJets_R04 = dataset("FullR04",NFIN=0,range=9,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',color=2,style=24,rebin=1)
  
  response, jetPt = MC_FullJets_R04.get2DHist('TrackJtCorrBin',dir='AliJJetJtTask/AliJJetJtMCHistManager', jetpt = True)
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  linear = ROOT.TF1("linear","x",0,10)
  linear.SetLineColor(1)
  linear.SetLineWidth(2)
  linear.SetLineStyle(2)
  x =linspace(0,10,50)
  y = linspace(0,10,50)
  

  #response[1].Draw('colz')
  #linear.Draw('Same')
  ij = 3
  rplt.hist2d(response[ij],label=MC_FullJets_R04.name(),norm=LogNorm(),colorbar=True)
  ax = plt.gca()
  ax.set_xlim([0,4])
  ax.set_ylim([0,4])
  ax.set_xlabel(r'$j_{T,true}\left[GeV\right]$') #Add x-axis labels for bottom row
  ax.set_ylabel(r'$j_{T,obs}\left[GeV\right]$') #Add x-axis labels for bottom row

  plt.plot(x,y,'-r')

  ax.text(0.5,2.8,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] +  '\n Jet Cone \n' + r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(jetPt[ij][0],jetPt[ij][1]),
          fontsize = 10,bbox=dict(boxstyle="round",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),
                   alpha = 0.8
                   ))
  
  plt.savefig("PythonFigures/ResponseMatrixNFin00JetPt03.pdf",format='pdf') #Save figure

  d['system'] = r'pPb MC$\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$'

  fig, axs = defs.makegrid(4,2,xlog=False,ylog=False,d=d,shareY=True,xtitle=r'$j_{T,true}\left[GeV\right]$',ytitle=r'$j_{T,obs}\left[GeV\right]$')
  axs = axs.reshape(8)
  axs[1].text(0.5,5,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] +  '\n Jet Cone',fontsize = 7,bbox=dict(boxstyle="round",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),
                   alpha=0.7
                   ))
  for r,pT,ax,i in zip(response[1:],jetPt[1:],axs,range(0,9)):
    rplt.hist2d(r,axes=ax,label=MC_FullJets_R04.name(),norm=LogNorm())
    x1 =linspace(0,10,50)
    y1 = linspace(0,10,50)
    plt.plot(x,y,'-r')
    
    ax.text(3.5,0.5,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1]),bbox=dict(boxstyle="round",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),
                   alpha=0.7
                   )
         )  
 
    ax.set_xlim([0,6.5]) #Set x-axis limits
    ax.set_ylim([0,6.5]) #Set y-axis limits
    
  axs[0].legend(loc = 'lower left')
       
  #plt.savefig("PythonFigures/MixedFullJetsR04JetConeJtSignal.pdf",format='pdf') #Save figure
  plt.savefig("PythonFigures/ResponseMatrixNFin00.pdf",format='pdf') #Save figure

  plt.show() #Draw figure on screen


  
if __name__ == "__main__": main()
