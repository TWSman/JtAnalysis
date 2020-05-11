import rootpy
import defs
import re
import matplotlib
from rootpy.plotting import Canvas
from rootpy.plotting import Legend
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
from dataset import *
def main(): 
  print('Number of arguments: ', len(sys.argv), 'arguments.')
  print('Argument list:',str(sys.argv))
  filenameOriginal = sys.argv[1]
  filenameNew = sys.argv[2]
  print("Input file Original: ")
  print(filenameOriginal)
  print("Input file new: ")
  print(filenameNew)
  #FullJets_R04 = dataset("FullR04",NFIN=0,filename=filename,directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=5)
  FullJets_R04_Original = dataset("Original",NFIN=0,range=8,filename=filenameOriginal,directory='AliJJetJtTask/AliJJetJtHistManager',color=2,style=24,rebin=5)
  FullJets_R04_New = dataset("New",NFIN = 0, range = 8, filename = filenameNew, directory = 'AliJJetJtTask/AliJJetJtHistManager',color=3, style=24,rebin = 5)
  ChargedJets_R04_Original = dataset("Original Charged",NFIN=2,range=8,filename=filenameOriginal,directory='AliJJetJtTask/AliJJetJtHistManager',color=2,style=24,rebin=5)
  ChargedJets_R04_New = dataset("New Charged",NFIN = 2, range = 8, filename = filenameNew, directory = 'AliJJetJtTask/AliJJetJtHistManager',color=3, style=24,rebin = 5)

  jetPt_Original = FullJets_R04_Original.getJetPt()
  jetPt_New = FullJets_R04_New.getJetPt()
  jetPt_Original.SetMarkerColor(1)
  jetPt_New.SetMarkerColor(2)
  
  ratio = jetPt_Original.Clone()
  ratio.Divide(jetPt_New)
  fig, axs = plt.subplots(2,1,figsize=(5,8),sharey=False,sharex=True)
  axs = axs.reshape(2)
  axs[0].text(15,1,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n Full jets\n' r'Anti-$k_T$, R=0.4',fontsize=7)
  axs[0].set_ylabel(r'$\frac{dN}{dp_{T}}$',fontsize=18)
  axs[0].set_xlabel(r'$p_{T}$',fontsize = 18)
  axs[1].set_xlabel(r'$p_{T}$',fontsize = 18)
  axs[1].set_ylabel('Ratio',fontsize=18)

  ax = axs[0]
  ax.set_xscale('log') #Set logarithmic scale
  ax.set_yscale('log')

  rplt.errorbar(jetPt_Original,xerr=False,emptybins=False,axes=ax,label='Original',fmt='+') #Plot jT histogram, 
  rplt.errorbar(jetPt_New,xerr=False,emptybins=False,axes=ax,label='New',fmt='go') #Plot bg jT histogram
  ax.legend(loc ='lower left')
  ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
  ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
  ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
  #ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
  ax.set_xlim([5,500]) #Set x-axis limits
  ax.set_ylim([5e-3,2e10]) #Set y-axis limits
  ax.grid(True) #Draw grid 

  ax= axs[1]
  rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,fmt='o')
  ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
  ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
  ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
  #ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
  ax.set_xscale('log') #Set logarithmic scale
  ax.set_yscale('log')
  ax.set_xlim([5,500]) #Set x-axis limits
  ax.set_ylim([1,10]) #Set y-axis limits
  ax.grid(True) #Draw grid  

  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
  plt.savefig('PythonFigures/JetPtComparison.pdf',format='pdf') #Save figure
  plt.show() #Draw figure on screen
  
if __name__ == "__main__": main()
