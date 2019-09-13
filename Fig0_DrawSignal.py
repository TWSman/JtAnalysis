import logging
from matplotlib.pyplot import box
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
import rootpy
import os.path
import defs
import re
import matplotlib
from rootpy.plotting import Canvas
from rootpy.plotting import Legend,Hist
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
from dataset import *
def main(): 

  start = 4
  end = 8
  n_figs = end-start
  title = "Full jets R=0.4"
  if(os.path.exists('RootFiles/Fig0.root')):
    inFile = "RootFiles/Fig0.root"
    inF = root_open(inFile,'r')
    signal = [inF.Get("jTSignalJetPt{:02d}".format(i)) for i in range(8)]
    jetPt = [(int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(1)),int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(2))) for h in signal] #Use regular expressions to extract jet pT range from histogram titles

  else:
    filename = "CF_pPb_legotrain/legotrain_CF_pPb_1839_20180613_LHC13bcde.root"

    print("Number of figs: {}".format(n_figs))
    print "Input file: "
    print filename
    Mixed_FullJets_R04 = datasetMixed(title,NFIN=0,range=(1,5),filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=2)
    signal,jetPt = Mixed_FullJets_R04.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = True)
    
    outFile = "Fig0.root"
    outF = root_open(outFile,"w+")
    for s,i in zip(signal,range(10)):
      s.SetName("jTSignalJetPt{:02d}".format(i))
      s.Write()
    outF.Close()
  
  
  
  n_rows = n_figs//4
  fig, axs = defs.makegrid(4,n_figs//4,xlog=True,ylog=True,d=d,shareY=True,figsize=(10,2.5))
  axs = axs.reshape(n_figs)
  axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone',fontsize = 7)
  for jT,pT,ax,i in zip(signal[start:],jetPt[start:],axs,range(0,9)):
    plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=title,fmt='o',fillstyle='none',ecolor='blue') #Plot jT histogram, 
    line = plot.get_children()[0]
    #line.set_markerfacecolor('none')
    #line.set_markeredgecolor('red')
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 

    ax.set_xlim([0.1,22]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
   
      
  plt.savefig("PythonFigures/MixedFullJetsR04JetConeJtSignalPtFrom{}To{}.pdf".format(start,end),format='pdf') #Save figure
  plt.show() #Draw figure on screen


  
if __name__ == "__main__": main()
