#!/usr/bin/python3
# dataset.py by Tomas Snellman
# Class to hold one dataset in Jet jT analysis
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
class dataset():
    
  def __init__(self,name,**kwargs):
    print("Creating {}".format(name))
    self._NFIN = kwargs['NFIN']
    self._filename = kwargs['filename']
    self._name = name
    self._directory = kwargs['directory']
    self.properties = kwargs
    self._f = root_open(self._filename)
    self._rebin = kwargs.get('rebin',1)
    self._measN = None
    self._measBgN = None
    self._range = kwargs.get('range',9)

    with root_open(self._filename,'read') as f:
      self._measN = [int(f.Get('{}/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}'.format(self._directory,self._NFIN,i)).GetEntries()) for i in range(0,self._range)] #Get number of jets by jet pT bins
      #print self._measN
      self._measBgN = [int(f.Get('{}/BgTrkNumberBin/BgTrkNumberBinNFin{:02d}JetPt{:02d}'.format(self._directory,self._NFIN,i)).GetEntries()) for i in range(0,self._range)] #Get number of background jets


  def printStats(self,**kwargs):
    ratios = [a / float(b) for b,a in zip(self._measN,self._measBgN)]
    hist = [self._f.Get('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory, 'histname':'JetPtBin','NFin':self._NFIN,'pT':i})) for i in range(0,self._range)]  #Get jT histograms from file an array
    jetPt = [(int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(1)),int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(2))) for h in hist] #Use regular expressions to extract jet pT range from histogram titles
    if(kwargs.get('format','') == 'latex'):
      if(kwargs.get('what','') == 'all'):
        print '\\begin{tabular}{lllllllllll}'
        print '{}'.format(self.name()),
      if(kwargs.get('what','') == 'all' or kwargs.get('what','') == 'jetpt' ):
        for pT in jetPt:
          print '& {}-{}'.format(pT[0],pT[1]),
        print '\\\\'
      if(kwargs.get('what','') == 'all' or kwargs.get('what','') == 'jets' ):
        if(kwargs.get('what','')=='all'):
          print 'Jets',
        else: print '{}'.format(self.name()),
        for N in self._measN:
            print '& {}'.format(N), 
        print '\\\\'
      if(kwargs.get('what','') == 'all' or kwargs.get('what','') == 'bg' ):
        if(kwargs.get('what','')=='all'):
          print 'Bg Cones',
        else: print '{}'.format(self.name()),
        for N in self._measBgN:
          print '& {}'.format(N), 
        print '\\\\'
      if(kwargs.get('what','') == 'all' or kwargs.get('what','') == 'bgratio' ):
        if(kwargs.get('what','')=='all'):
          print 'Bg ratio',
        else: print '{}'.format(self.name()),
        for N in ratios:
          print '& {:.2f}\\%'.format(100*N), 
        print '\\\\'
      if(kwargs.get('what','') == 'all'):
        print '\\end{tabular}'
    else:
      print "{} Jets: ".format(self.name())
      print self._measN
      print "{} Bg Cones: ".format(self.name())
      print self._measBgN
      print "{} Bg Ratio: ".format(self.name())
      print ratios

        
     
  def getHist(self,name,**kwargs):
    hist = [self._f.Get('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory, 'histname':name,'NFin':self._NFIN,'pT':i})) for i in range(0,self._range)]  #Get jT histograms from file an array
    #print('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory, 'histname':name,'NFin':self._NFIN,'pT':1}))
    jetPt = [(int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(1)),int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(2))) for h in hist] #Use regular expressions to extract jet pT range from histogram titles
    #print(len(hist))
    #print hist
    #print jetPt
    for h,N,bgN in zip(hist,self._measN,self._measBgN):
      h.Sumw2()
      print "Rebinning {} by {} in set {} that has {} bins".format(h.GetTitle(),self._rebin,self._name,h.GetNbinsX())
      h.Rebin(self._rebin)
      if(kwargs.get('isBg',False)):
        h.SetLineColor(self.properties.get('color',1) + 1)
        h.SetMarkerColor(self.properties.get('color',1) + 1)
        h.Scale(1.0/bgN,'width')
        print("{} is bg".format(name))

      else:
        h.SetLineColor(self.properties.get('color',1))
        h.SetMarkerColor(self.properties.get('color',1))
        h.Scale(1.0/N,'width')

      h.SetMarkerStyle(self.properties.get('style',24))
      h.SetMarkerSize(0.5)
      h.SetLineColor(1)

    if(kwargs.get('jetpt',False)):
      return hist,jetPt
    else:
      return hist
  
  def name(self):
    return self._name
  
  def getprop(self,key):
    return self.properties.get(key,None)

d = dict(
    jetalg = r'Anti-$k_T$, R=0.4',
    jettype = 'Charged Jets',
    system = r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' ,
    trigger = 'kEMCEJE'
)

def compareSets(sets,histname):
  datasets = [dataset.getHist(histname,jetpt = True) for dataset in sets]
  names = [dataset.name() for dataset in sets]
  fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d)
  axs = axs.reshape(8)
  axs[1].text(0.2,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'],fontsize = 7)
  
  for i,(set,name) in enumerate(zip(datasets,names)):
    for jT,pT,ax,j in zip(set[0][1:],set[1][1:],axs,range(0,9)):
      rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=name,fmt='o') #Plot jT histogram, 
      if(i == 0):
        ax.set_xlim([0.1,20]) #Set x-axis limits
        ax.set_ylim([5e-4,2e3]) #Set y-axis limits
        
  axs[0].legend(loc = 'lower left')
  plt.show() #Draw figure on screen

def compareSetsWithRatio(sets,histname):
  datasets = [dataset.getHist(histname,jetpt = True) for dataset in sets]
  names = [dataset.name() for dataset in sets]
  fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d)
  axs = axs.reshape(8)
  axs[1].text(0.2,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'],fontsize = 7)
  
  for i,(set,name) in enumerate(zip(datasets,names)):
    for jT,pT,ax,j in zip(set[0][1::2],set[1][1::2],axs[0:4],range(0,5)):
      rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=name,fmt='o') #Plot jT histogram, 
      ax.set_xlim([0.01,20]) #Set x-axis limits
      ax.set_ylim([5e-4,2e3]) #Set y-axis limits
      ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 

  
  ratiosets = []
  
  for set in datasets:
    ratios = []
    for s,divider in zip(set[0],datasets[0][0]):
      h = s.Clone()
      h.Divide(divider)
      ratios.append(h)
    ratios1 = (ratios,set[1])
    ratiosets.append(ratios1)
    
  for i,(set,name) in enumerate(zip(ratiosets[1:],names[1:])):
    for ratio,pT,ax,j in zip(set[0][1::2],set[1][1::2],axs[4:9],range(0,5)):
      rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label=name,fmt='o') #Plot ratio histogram, 
      #if(i == 0):
      ax.set_xlim([0.01,20]) #Set x-axis limits
      ax.set_ylim([0.1,2]) #Set y-axis limits     
      ax.set_yscale('linear')
  axs[0].legend(loc = 'lower left')
  plt.show() #Draw figure on screen
  
def JtWithBackgroundRatio(dataset,histname,bgname,all=False):
  jtHistos = dataset.getHist(histname,jetpt = True) 
  bgHistos = dataset.getHist(bgname,jetpt = True,isBg = True)
  fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d)
  axs = axs.reshape(8)
  axs[1].text(0.02,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'] + '\n' + dataset.name(),fontsize = 7)
  
  for jT,bgJt,pT,ax,j in zip(jtHistos[0][1::2],bgHistos[0][1::2],jtHistos[1][1::2],axs[0:4],range(0,5)):
    rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label='Jt',fmt='o') #Plot jT histogram,
    rplt.errorbar(bgJt,xerr=False,emptybins=False,axes=ax,label='BgJt',fmt='o') #Plot jT histogram, 
    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 

  
  ratios = []
  for jT,bg in zip(jtHistos[0],bgHistos[0]):
    h = bg.Clone()
    h.Divide(jT)
    ratios.append(h)

  axs[4].set_ylabel('Ratio \n Inclusive/Background',fontsize=9) #Add y-axis labels to left- and righmost subfigures
  axs[-1].set_ylabel('Ratio \n Inclusive/Background',fontsize=9)
  for ratio,pT,ax,j in zip(ratios[1::2],jtHistos[1][1::2],axs[4:9],range(0,5)):
    rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram, 
    #if(i == 0):
    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([0,2.2]) #Set y-axis limits     
    ax.set_yscale('linear')
  axs[0].legend(loc = 'lower left')
  plt.show() #Draw figure on screen
  
def JtWithBackgroundRatioAll(dataset,histname,bgname,all=False):
  jtHistos = dataset.getHist(histname,jetpt = True) 
  bgHistos = dataset.getHist(bgname,jetpt = True,isBg = True)
  fig, axs = defs.makegrid(4,4,xlog=True,ylog=True,d=d)
  axs = axs.reshape(16)
  axs[1].text(0.02,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'] + '\n' + dataset.name(),fontsize = 7)
  print(type(axs))
  for jT,bgJt,pT,ax,j in zip(jtHistos[0][1:],bgHistos[0][1:],jtHistos[1][1:],axs[0:4].tolist()+axs[8:12].tolist(),range(0,8)):
    rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label='Jt',fmt='o') #Plot jT histogram,
    rplt.errorbar(bgJt,xerr=False,emptybins=False,axes=ax,label='BgJt',fmt='o') #Plot jT histogram, 
    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 

  
  ratios = []
  for jT,bg in zip(jtHistos[0],bgHistos[0]):
    h = bg.Clone()
    h.Divide(jT)
    ratios.append(h)

  axs[4].set_ylabel('Ratio \n Inclusive/Background',fontsize=9) #Add y-axis labels to left- and righmost subfigures
  axs[-1].set_ylabel('Ratio \n Inclusive/Background',fontsize=9)
  for ratio,pT,ax,j in zip(ratios[1:],jtHistos[1][1:],axs[4:8].tolist()+axs[12:16].tolist(),range(0,8)):
    rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram, 
    #if(i == 0):
    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([0,2.2]) #Set y-axis limits     
    ax.set_yscale('linear')
  axs[0].legend(loc = 'lower left')
  plt.show() #Draw figure on screen 
    
def main():
  print 'Number of arguments: ', len(sys.argv), 'arguments.'
  print 'Argument list:',str(sys.argv)
  filename = sys.argv[1]
  print "Input file: "
  print filename
  MB_FullJets_R04 = dataset("MBFullR04",NFIN=0,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',color=1,style=24,rebin=5)
  MB_FullJets_R05 = dataset("MBFullR05",NFIN=1,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',color=2,style=24,rebin=5)
  MB_ChargedJets_R04 = dataset("MBChargedR04",NFIN=2,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',color=1,style=24,rebin=5)
  MB_ChargedJets_R05 = dataset("MBChargedR05",NFIN=3,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',color=2,style=24,rebin=5,range = 8)
  triggered_FullJets_R04 = dataset("TriggeredFullR04",NFIN=0,filename=filename,directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=5)
  triggered_FullJets_R05 = dataset("TriggeredFullR05",NFIN=1,filename=filename,directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=1,style=24,rebin=5)
  triggered_ChargedJets_R04 = dataset("TriggeredChargedR04",NFIN=2,filename=filename,directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=3,style=24,rebin=5)
  triggered_ChargedJets_R05 = dataset("TriggeredChargedR05",NFIN=3,filename=filename,directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=4,style=24,rebin=5)
  datasets = (MB_FullJets_R04,MB_FullJets_R05,MB_ChargedJets_R04,MB_ChargedJets_R05,triggered_FullJets_R04,triggered_FullJets_R05,triggered_ChargedJets_R04,triggered_ChargedJets_R05)
  for set in datasets:
    set.printStats(format ='latex',what = 'jets')
  for set in datasets:
    set.printStats(format ='latex',what = 'bg')
  for set in datasets:
    set.printStats(format ='latex',what = 'bgratio')
 

  #measJt,jetpt = MB_FullJets_R04.getHist('JetConeJtWeightBin',jetpt = True)
  
  #compareSetsWithRatio((MB_FullJets_R05,MB_FullJets_R04),'JetConeJtWeightBin')
  #compareSetsWithRatio((triggered_FullJets_R05,triggered_FullJets_R04),'JetConeJtWeightBin')
  #compareSetsWithRatio((triggered_ChargedJets_R04,triggered_FullJets_R04),'JetConeJtWeightBin')
  JtWithBackgroundRatioAll(MB_FullJets_R04, 'JetConeJtWeightBin', 'BgJtWeightBin')


#===============================================================================
#   fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d)
#   axs = axs.reshape(8)
#   axs[1].text(0.2,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'],fontsize = 7)
# 
#   for jT,pT,ax,i in zip (measJt[1:],jetpt[1:],axs,range(0,9)):
#     rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label='jT',fmt='+') #Plot jT histogram, 
#     ax.set_xlim([0.1,20]) #Set x-axis limits
#     ax.set_ylim([5e-4,2e3]) #Set y-axis limits
#     #if i == 0: #For the first subfigure add a legend to bottom left corner
#     #  ax.legend(loc ='lower left')
#   #plt.savefig("jTJetPtBins",format='pdf') #Save figure
#   axs[0].legend(loc = 'lower left')
#   plt.show() #Draw figure on screen
#===============================================================================

  

    


if __name__ == "__main__": main()
    
        