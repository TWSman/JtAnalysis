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
  print("Input file: ")
  print(filename)
  #FullJets_R04 = dataset("FullR04",NFIN=0,filename=filename,directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=5)
  Mixed_FullJets_R04 = datasetMixed("FullR04",NFIN=0,range=5,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=5)
  
  bgname = 'BgJtWeightBin';
  rndmbgname = 'BgRndmJtWeightBin';
  print("Get Inclusive histograms")
  inclusive,jetPt = Mixed_FullJets_R04.getHist('JetConeJtWeightBin',jetpt = True)
  print("Get random background histograms")
  rndmbgHistos = Mixed_FullJets_R04.getHist(rndmbgname,jetpt = True,isRndmBg = True) 
  print("Get perpendicular cone background histograms")
  perpbgHistos = Mixed_FullJets_R04.getHist(bgname,jetpt = True,isBg = True)
  
  
  print("Get perpendicular cone bg subtracted signal histograms")
  signalPerp,jetPt = Mixed_FullJets_R04.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = True)
  print("Get random bg subtracted signal histograms")

  signalRand,jetPt = Mixed_FullJets_R04.getSubtracted('JetConeJtWeightBin','BgRndmJtWeightBin',jetpt = True,randomBG=True)  


  fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d,shareY = False)
  axs = axs.reshape(8)
  axs[1].text(0.02,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'] + '\n' + Mixed_FullJets_R04.name(),fontsize = 7)
  
  for rand,perp,pT,ax,j in zip(rndmbgHistos[0][1::2],perpbgHistos[0][1::2],rndmbgHistos[1][1::2],axs[0:4],range(0,5)):
    rplt.errorbar(rand,xerr=False,emptybins=False,axes=ax,label='Random BG',fmt='o') #Plot rand histogram,
    rplt.errorbar(perp,xerr=False,emptybins=False,axes=ax,label='Perp. Cone BG',fmt='o') #Plot perp histogram, 
    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 

  
  ratios = []
  for rand,bg in zip(rndmbgHistos[0],perpbgHistos[0]):
    h = bg.Clone()
    h.Divide(rand)
    ratios.append(h)

  axs[4].set_ylabel('Ratio \n Perp/Random',fontsize=9) #Add y-axis labels to left- and righmost subfigures
  axs[-1].set_ylabel('Ratio \n Perp/Random',fontsize=9)
  for ratio,pT,ax,j in zip(ratios[1::2],rndmbgHistos[1][1::2],axs[4:9],range(0,5)):
    rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram, 
    #if(i == 0):
    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([0.43,1.3]) #Set y-axis limits     
    ax.set_yscale('linear')
  axs[0].legend(loc = 'lower left')

  plt.savefig("PythonFigures/MixedFullJetsR04BackgroundComparison.pdf",format='pdf') #Save figure
  plt.show() #Draw figure on screen
  
  fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d,shareY = False)
  axs = axs.reshape(8)
  axs[1].text(0.02,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'] + '\n' + Mixed_FullJets_R04.name(),fontsize = 7)
  
  for rand,perp,randBg,inc,pT,ax,j in zip(signalRand[1::2],signalPerp[1::2],rndmbgHistos[0][1::2],inclusive[1::2],jetPt[1::2],axs[0:4],range(0,5)):
    rand.SetMarkerColor(2)
    perp.SetMarkerColor(3)
    randBg.SetMarkerColor(1)
    inc.SetMarkerColor(4)
    rplt.errorbar(rand,xerr=False,emptybins=False,axes=ax,label='Random BG signal',fmt='o') #Plot rand histogram,
    rplt.errorbar(perp,xerr=False,emptybins=False,axes=ax,label='Perp. Cone BG signal',fmt='o') #Plot jT histogram, 



    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 

  
  ratios = []
  for rand,perp in zip(signalRand,signalPerp):
    h = perp.Clone()
    h.Divide(rand)
    ratios.append(h)

  axs[4].set_ylabel('Ratio \n Perp/Random',fontsize=9) #Add y-axis labels to left- and righmost subfigures
  axs[-1].set_ylabel('Ratio \n Perp/Random',fontsize=9)
  for ratio,pT,ax,j in zip(ratios[1::2],jetPt[1::2],axs[4:9],range(0,5)):
    rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram, 
    #if(i == 0):
    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([0.85,1.3]) #Set y-axis limits     
    ax.set_yscale('linear')
  axs[0].legend(loc = 'lower left')

  plt.savefig("PythonFigures/MixedFullJetsR04SignalBackgroundComparison.pdf",format='pdf') #Save figure
  plt.show() #Draw figure on screen
  
  fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d,shareY = False)
  axs = axs.reshape(8)
  axs[1].text(0.02,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'] + '\n' + Mixed_FullJets_R04.name(),fontsize = 7)
  
  for rand,perp,randBg,inc,pT,ax,j in zip(signalRand[1::2],signalPerp[1::2],rndmbgHistos[0][1::2],inclusive[1::2],jetPt[1::2],axs[0:4],range(0,5)):
    rand.SetMarkerColor(2)
    perp.SetMarkerColor(3)
    randBg.SetMarkerColor(1)
    inc.SetMarkerColor(4)
    rplt.errorbar(rand,xerr=False,emptybins=False,axes=ax,label='Random BG signal',fmt='o') #Plot rand histogram,
    rplt.errorbar(perp,xerr=False,emptybins=False,axes=ax,label='Perp. Cone BG signal',fmt='o') #Plot jT histogram, 
    rplt.errorbar(randBg,xerr=False,emptybins=False,axes=ax,label='Random BG',fmt='o') #Plot jT histogram, 
    rplt.errorbar(inc,xerr=False,emptybins=False,axes=ax,label='Inclusive',fmt='o') #Plot jT histogram, 


    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 

  
  ratiosRand = []
  ratiosPerp = []
  ratiosBg = []
  for rand,perp,bg,inc in zip(signalRand,signalPerp,rndmbgHistos[0],inclusive):
    h1 = rand.Clone()
    h1.Divide(inc)
    ratiosRand.append(h1)
    
    h2 = perp.Clone()
    h2.Divide(inc)
    ratiosPerp.append(h2)
    
    h3 = bg.Clone()
    h3.Divide(inc)
    ratiosBg.append(h3)

  axs[4].set_ylabel('Ratio \n Perp/Random',fontsize=9) #Add y-axis labels to left- and righmost subfigures
  axs[-1].set_ylabel('Ratio \n Perp/Random',fontsize=9)
  for ratioRand,ratioPerp,ratioBg,pT,ax,j in zip(ratiosRand[1::2],ratiosPerp[1::2],ratiosBg[1::2],jetPt[1::2],axs[4:9],range(0,5)):
    ratioRand.SetMarkerColor(2)
    ratioPerp.SetMarkerColor(3)
    ratioBg.SetMarkerColor(1)
    rplt.errorbar(ratioRand,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram,
    rplt.errorbar(ratioBg,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram, 
    rplt.errorbar(ratioBg,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram, 
 
    #if(i == 0):
    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([0,1.05]) #Set y-axis limits     
    ax.set_yscale('linear')
  axs[0].legend(loc = 'lower left')

  plt.savefig("PythonFigures/MixedFullJetsR04QABackgroundComparison.pdf",format='pdf') #Save figure
  plt.show() #Draw figure on screen



  
if __name__ == "__main__": main()
