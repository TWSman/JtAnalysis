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
  separate = int(sys.argv[2])
  if(len(sys.argv) > 3):
    start = int(sys.argv[3])
    end = int(sys.argv[4])
  else:
    start = 1
    end = 9
  n_figs = end-start
  print("Number of figs: {}".format(n_figs))
  print "Input file: "
  print filename
  Mixed_FullJets_R04 = datasetMixed("Full jets R=0.4",NFIN=0,range=(1,5),filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=2)
  signal,jetPt = Mixed_FullJets_R04.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = True)
  if(False):
    if(separate > 0):
      fig = plt.figure(1)
      ax = fig.add_subplot(1,1,1)
      ax.set_xscale('log')
      ax.set_yscale('log')
      ax.text(0.2,0.0005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone',fontsize = 7)
      rplt.errorbar(signal[separate],xerr=False,emptybins=False,axes=ax,label=Mixed_FullJets_R04.name(),fmt='o') #Plot jT histogram, 
      ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(jetPt[separate][0],jetPt[separate][1])) 
      ax.set_xlim([0.1,12])
      ax.set_ylim([5e-6,1.5e3])
      ax.legend(loc = 'lower left')
      
      plt.savefig("PythonFigures/MixedFullJetsR04JetConeJtSignalJetPt{0}.pdf".format(separate),format='pdf') #Save figure
      plt.show() #Draw figure on screen
  
    else:
      n_rows = n_figs//4
      print(n_rows)
      fig, axs = defs.makegrid(4,n_figs//4,xlog=True,ylog=True,d=d,shareY=True,figsize=(10,2.5))
      axs = axs.reshape(n_figs)
      axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone',fontsize = 7)
      for jT,pT,ax,i in zip(signal[start:],jetPt[start:],axs,range(0,9)):
        plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=Mixed_FullJets_R04.name(),fmt='o',fillstyle='none',ecolor='blue') #Plot jT histogram, 
        line = plot.get_children()[0]
        #line.set_markerfacecolor('none')
        #line.set_markeredgecolor('red')
        ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
    
        ax.set_xlim([0.1,22]) #Set x-axis limits
        ax.set_ylim([5e-4,2e3]) #Set y-axis limits
        ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
       
      #axs[0].legend(loc = 'lower left')
          
      plt.savefig("PythonFigures/MixedFullJetsR04JetConeJtSignalPtFrom{}To{}.pdf".format(start,end),format='pdf') #Save figure
      plt.show() #Draw figure on screen

  fig,axs = defs.makeRatio(xlog=True,ylog=True,d=d,shareY=False,figsize=(10,10))
  axs = axs.reshape(2)
  ax = axs[0]
  for jT,pT,i in zip(signal[start:],jetPt[start:],range(0,9)):
    h = jT.Clone()
    h.Scale(10**i)
    h.SetMarkerColor(i)
    label = r'${:02d}\:\mathrm{{GeV}} < p_{{T,\mathrm{{jet}}}} < {:02d}\:\mathrm{{GeV}} \left(\times 10^{:01d}\right)$'.format(pT[0],pT[1],i)
    plot = rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,label=label,fmt='o',fillstyle='none') #Plot jT histogram,
  
  ax.set_xlim([0.1,22])
  ax.set_ylim([5e-4,2e8])
  ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
  ax.legend(loc = 'lower left')

  plt.show() #Draw figure on screen

  
if __name__ == "__main__": main()
