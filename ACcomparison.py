import logging
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
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
  filename = sys.argv[1]
  comment = sys.argv[2]
  if(len(sys.argv) > 3):
    start = int(sys.argv[3])
  else:
    start = 1
  print("Input file: ")
  print(filename)
  if comment in ["LHC13d","LHC13e","LHC13EMCAL"]:
    taskbase = "AliJJetJtTask_kEMCEJE"
  else:
    taskbase = "AliJJetJtTask"
  Aside_FullJets = dataset("Full jets R=0.4, A side",NFIN=0,range=(0,8),filename=filename,directory='{}_Aside/AliJJetJtHistManager'.format(taskbase),color=2,style=24,rebin=2)
  Cside_FullJets = dataset("Full jets R=0.4, C side",NFIN=0,range=(0,8),filename=filename,directory='{}_Cside/AliJJetJtHistManager'.format(taskbase),color=2,style=24,rebin=2)

  signal,jetPt = Aside_FullJets.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = True)
  signal2 = Cside_FullJets.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = False)
  signals = [signal,signal2]
  names = ["A side","C side"]
  colors = [1,2]
  n_figs = 8
  n_rows = 2
  fig, axs = defs.makegrid(4,n_figs//4,xlog=True,ylog=True,d=d,shareY=False,figsize=(15,7.5))
  axs = axs.reshape(8)
  axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone' + '\n Inclusive jT' + '\n' + comment,fontsize = 7)
  for signal,name,color,j in zip(signals,names,colors,range(10)):
    print("Plot {}".format(name))
    for jT,pT,ax,i in zip(signal[start:],jetPt[start:],axs[0:4],range(0,9)):
      print(i)
      jT.SetMarkerColor(color)
      jT.SetLineColor(color)
      plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=name,fmt='o',fillstyle='none',ecolor='blue') #Plot jT histogram, 
      ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
  
      ax.set_xlim([0.1,22]) #Set x-axis limits
      ax.set_ylim([5e-4,2e3]) #Set y-axis limits
      ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
    ratios = []
    if j > 0:
      for jT,div in zip(signal,signals[0]):
        h = jT.Clone()
        h.Divide(div)
        ratios.append(h)
      axs[4].set_ylabel('Ratio C/A',fontsize=9) #Add y-axis labels to left- and righmost subfigures
      axs[-1].set_ylabel('Ratio C/A',fontsize=9)       
      for ratio,pT,ax in zip(ratios[start:],jetPt[start:],axs[4:9]):
        rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram,
        print("Draw {}".format(ratio.GetName()))   
        #if(i == 0):
        ax.set_yscale('linear')
        ax.set_xlim([0.1,20]) #Set x-axis limits
        ax.set_ylim([0.8,1.2]) #Set y-axis limits
     
    axs[0].legend(loc = 'lower left')
        
  plt.savefig("PythonFigures/ACsideComparison/ACsideJetConeJtInclusivePtFrom{}To{}{}.pdf".format(start,start+4,comment),format='pdf') #Save figure
  plt.show() #Draw figure on screen


  
if __name__ == "__main__": main()
