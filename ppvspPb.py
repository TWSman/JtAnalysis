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
  filename1 = '~/OneDrive/work/032.JTAnalysis/Unfolding/RooUnfold/CF_pPb_legotrain/legotrain_CF_pPb_CF_pPb-1209_20170807_LHC13bcde.root'
  filename2 = '~/OneDrive/work/032.JTAnalysis/Unfolding/RooUnfold/CF_pp_legotrain/legotrain_CF_pp_1405_20170810-7TeV_LHC10_p2_AOD147.root'

  print "Input file: "
  print filename1
  Mixed_FullJets_R04 = datasetMixed("FullR04",NFIN=0,range=5,filename=filename1,directory='AliJJetJtTask/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=5)
  pp_FullJets_R04 = dataset("ppFullR04",NFIN=0,filename=filename2,directory='AliJJetJtTask/AliJJetJtHistManager',color=1,style=24,rebin=5)
  fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d,shareY=True)
  signal,jetPt = Mixed_FullJets_R04.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = True)
  axs = axs.reshape(8)
  axs[1].text(0.02,0.0005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'] + '\n' + Mixed_FullJets_R04.name() + '\n Jet Cone',fontsize = 7)
  for jT,pT,ax,i in zip(signal[1:],jetPt[1:],axs,range(0,9)):
    rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=Mixed_FullJets_R04.name(),fmt='o') #Plot jT histogram, 

    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 

    if(i == 0):
      ax.set_xlim([0.01,20]) #Set x-axis limits
      ax.set_ylim([5e-4,2e3]) #Set y-axis limits
   
  axs[0].legend(loc = 'lower left')
      
  #plt.savefig("PythonFigures/MixedFullJetsR04JetConeJtSignal.pdf",format='pdf') #Save figure
  plt.show() #Draw figure on screen


  
if __name__ == "__main__": main()
