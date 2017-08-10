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
  print 'Number of arguments: ', len(sys.argv), 'arguments.'
  print 'Argument list:',str(sys.argv)
  filename = sys.argv[1]
  print "Input file: "
  print filename
  Mixed_FullJets_R04 = datasetMixed("FullR04",NFIN=0,range=5,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=5)
  fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d,shareY=True)

  jTHistos, jetPt = Mixed_FullJets_R04.getHist('JetConeJtWeightBin',jetpt = True)
  bgHistos = Mixed_FullJets_R04.getHist('BgJtWeightBin', jetpt = False)
  signal = [h.Clone() for h in jTHistos]
  for (jT,Bg) in zip (signal,bgHistos):
    jT.Add(Bg,-1)
  
  axs = axs.reshape(8)
  axs[1].text(0.02,0.0005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'] + '\n' + Mixed_FullJets_R04.name() + '\n Jet Cone',fontsize = 7)
  for jT,pT,ax,i in zip(signal[1:],jetPt[1:],axs,range(0,9)):
    rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=Mixed_FullJets_R04.name(),fmt='o') #Plot jT histogram, 
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 

    if(i == 0):
      ax.set_xlim([0.01,20]) #Set x-axis limits
      ax.set_ylim([5e-4,2e3]) #Set y-axis limits
   
  axs[0].legend(loc = 'lower left')
      
  plt.savefig("PythonFigures/MixedFullJetsR04JetConeJtSignal.pdf",format='pdf') #Save figure
  plt.show() #Draw figure on screen


  
if __name__ == "__main__": main()
