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

class dataset(object):
  """ Basic class for handling a single dataset with one input file/directory and one jet finder
  
  Can retrieve histograms binned by jet Pt.
  
  Author Tomas Snellman tomas.snellman@cern.ch 
    
  Attributes:
    properties: List of kwargs used to create class
    _filename: Name of used file
    _NFIN: Index of jet finder
    _name: name of class. Used in figure legends
    _directory: Directory in the file where histograms are located
  """
  def __init__(self,name,**kwargs):
    """
    The constructor.
  
    Args:
      name: The name of the dataset, will show in figure legends  
      
    Kwargs:
      **filename** Name of the file\n
      **directory** Directory where the histograms are to be found\n
      **NFIN** Jet finder index\n
      **rebin** How much rebinning should be done\n
      **range** Which histograms should be retrieved, 9 by default which gives the full jet pT range from 5-10GeV bin to 150-150GeV bin\n 
    """
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
      self._measBgN = [int(f.Get('{}/BgTrkNumberBin/BgTrkNumberBinNFin{:02d}JetPt{:02d}'.format(self._directory,self._NFIN,i)).GetEntries()) for i in range(0,self._range)] #Get number of background jets

  def printStats(self,**kwargs):
    """Print available statistics, number of jets by jet PT bin
    
    Args:
    
    Kwargs:
    **format** if set to latex will print a latex table source code, otherwise python lists \n
    **what** What statistics to give \n
    jets: gives number of jets\n
    bg: gives number of background cones \n
    bgratio: Gives ratio of bg cones to number of jets \n
    """   
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
   

  def getSubtracted(self,inclusive,background,**kwargs):
    hist,jetpt = self.getHist(inclusive,jetpt=True)
    bg = self.getHist(background,isBG = True,jetpt = False)
    signals = []
    for h,b in zip(hist,bg):
      s = h.Clone()
      s.Add(b,-1)
      signals.append(s)
    if(kwargs.get('jetpt',False)):
      return signals,jetpt
    else:
      return signals    
  
  def getHist(self,name,**kwargs):
    """
    Retrieve a list of histograms by jet pT bins
    
    Args:
      name: Name of histogram
  
    Kwargs:
      isbg: Determines which normalization to use    
    
    """
    if('dir' in kwargs):
      hist = [self._f.Get('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':kwargs['dir'], 'histname':name,'NFin':self._NFIN,'pT':i})).Clone() for i in range(0,self._range)]  #Get jT histograms from file an array
    else:
      hist = [self._f.Get('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory, 'histname':name,'NFin':self._NFIN,'pT':i})).Clone() for i in range(0,self._range)]  #Get jT histograms from file an array
    #print('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory, 'histname':name,'NFin':self._NFIN,'pT':1}))
    jetPt = [(int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(1)),int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(2))) for h in hist] #Use regular expressions to extract jet pT range from histogram titles
    #print(len(hist))
    #print hist
    #print jetPt
    for h,N,bgN in zip(hist,self._measN,self._measBgN):
      h.Sumw2()
      print "Rebinning {} by {} in set {} that has {} bins".format(h.GetTitle(),self._rebin,self._name,h.GetNbinsX())
      h.Rebin(self._rebin)
      print(kwargs)
      if(self.properties.get('isWeight',False)):
        h.SetLineColor(self.properties.get('color',1))
        h.SetMarkerColor(self.properties.get('color',1))
        h.Scale(1.0,'width')
      else:
        if(kwargs.get('isBg', False)):
          h.SetLineColor(self.properties.get('color', 1) + 1)
          h.SetMarkerColor(self.properties.get('color', 1) + 1)
          h.Scale(1.0 / bgN, 'width')
          print("{} is bg".format(name))
        else:
          h.SetLineColor(self.properties.get('color', 1))
          h.SetMarkerColor(self.properties.get('color', 1))
          h.Scale(1.0 / N, 'width')

      h.SetMarkerStyle(self.properties.get('style',24))
      h.SetMarkerSize(0.5)
      h.SetLineColor(1)

    if(kwargs.get('jetpt',False)):
      return hist,jetPt
    else:
      return hist
    
  def get2DHist(self,name,**kwargs):
    """
    Retrieve a list of 2D histograms by jet pT bins
    
    Args:
      name: Name of histogram
  
    Kwargs:
      isbg: Determines which normalization to use    
    
    """
    if('dir' in kwargs):
      hist = [self._f.Get('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':kwargs['dir'], 'histname':name,'NFin':self._NFIN,'pT':i})).Clone() for i in range(0,self._range)]  #Get jT histograms from file an array
    else:
      hist = [self._f.Get('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory, 'histname':name,'NFin':self._NFIN,'pT':i})).Clone() for i in range(0,self._range)]  #Get jT histograms from file an array
    jetPt = [(int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(1)),int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(2))) for h in hist] #Use regular expressions to extract jet pT range from histogram titles
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
    jettype = 'Full Jets',
    system = r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' ,
    trigger = 'kEMCEJE'
)

class datasetMixed(dataset,object):
  """Class for handling datasets where different bins come from different files. For example MB data for low jet pT and triggered data for high jet pT
  
  """
  def __init__(self,name,**kwargs):
    """Constructor
    
    """
    super(datasetMixed,self).__init__(name,**kwargs)
    self._filename1 = kwargs['filename']
    self._filename2 = kwargs.get('filename2',self._filename1)
    if(self._filename2 != self._filename1):
      self._f2 = root_open(self._filename)
    else: 
      self._f2 = self._f
    self._directory1 = kwargs['directory']
    self._directory2 = kwargs['directory2']
    print "Directory 1: {} Directory 2: {}".format(self._directory1,self._directory2)
    self._measN1 = [int(self._f.Get('{}/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}'.format(self._directory1,self._NFIN,i)).GetEntries()) for i in range(0,self._range)] #Get number of jets by jet pT bins
    self._measN2 = [int(self._f2.Get('{}/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}'.format(self._directory2,self._NFIN,i)).GetEntries()) for i in range(0,9)] #Get number of jets by jet pT bins
    self._measBgN1 = [int(self._f.Get('{}/BgTrkNumberBin/BgTrkNumberBinNFin{:02d}JetPt{:02d}'.format(self._directory1,self._NFIN,i)).GetEntries()) for i in range(0,self._range)] #Get number of background jets
    self._measBgN2 = [int(self._f2.Get('{}/BgTrkNumberBin/BgTrkNumberBinNFin{:02d}JetPt{:02d}'.format(self._directory2,self._NFIN,i)).GetEntries()) for i in range(0,9)] #Get number of background jets
    self._measN = [self._measN1[i] if (i < self._range) else self._measN2[i] for i in range(0,9)]
    self._measBgN = [self._measBgN1[i] if (i < self._range) else self._measBgN2[i] for i in range(0,9)]
    print "MeasN1: {} MeasN2: {}".format(self._measN1,self._measN2)
    print "MeasN: {} ".format(self._measN)
  
  def getHist(self,name,**kwargs):
    """Return a set of histograms
    
    Args:
      self: Pointer to the object
      name: Name of histograms
      kwargs: list of keyword arguments
    """
    hist1 = [self._f.Get('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory1, 'histname':name,'NFin':self._NFIN,'pT':i})).Clone() for i in range(0,self._range)]  #Get jT histograms from file an array
    hist2 = [self._f2.Get('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory2, 'histname':name,'NFin':self._NFIN,'pT':i})).Clone() for i in range(0,9)]  #Get jT histograms from file an array
    hist = [hist1[i] if(i < self._range) else hist2[i] for i in range(0,9)]    
    #print('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory, 'histname':name,'NFin':self._NFIN,'pT':1}))
    jetPt = [(int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(1)),int(re.search( r'p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*',h.GetTitle(), re.M|re.I).group(2))) for h in hist] #Use regular expressions to extract jet pT range from histogram titles
    #print(len(hist))
    #print hist
    #print jetPt
    for h,N,bgN in zip(hist,self._measN,self._measBgN):
      h.Sumw2()
      print "Rebinning {} by {} in set {} that has {} bins".format(h.GetTitle(),self._rebin,self._name,h.GetNbinsX())
      h.Rebin(self._rebin)
      if(self.properties.get('isWeight',False)):
        h.SetLineColor(self.properties.get('color',1))
        h.SetMarkerColor(self.properties.get('color',1))
        h.Scale(1.0,'width')
      else:
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
    
def compareSets(sets,histname):
  """Draws a comparison between sets of a single histogram
  
  Args:
    sets: List of sets to be used
    histname: Name of histogram to be drawn
  """
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
  """Draw a comparison between sets in 4 jet pT bins with ratios to first set in list
  
  Args:
    sets: List os sets to be used
    histname: Name of histogram to be plotted
    
  """
  datasets = [dataset.getHist(histname,jetpt = True) for dataset in sets]
  names = [dataset.name() for dataset in sets]
  fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d,shareY = False)
  axs = axs.reshape(8)
  axs[1].text(0.02,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'],fontsize = 7)
  
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

  axs[4].set_ylabel('Ratio',fontsize=9) #Add y-axis labels to left- and righmost subfigures
  axs[-1].set_ylabel('Ratio',fontsize=9)
  
  for i,(set,name) in enumerate(zip(ratiosets[1:],names[1:])):
    for ratio,pT,ax,j in zip(set[0][1::2],set[1][1::2],axs[4:9],range(0,5)):
      rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label=name,fmt='o') #Plot ratio histogram, 
      #if(i == 0):
      ax.set_xlim([0.01,20]) #Set x-axis limits
      ax.set_ylim([0.1,1.8]) #Set y-axis limits     
      ax.set_yscale('linear')
  axs[0].legend(loc = 'lower left')
  
  
def JtWithBackgroundRatio(dataset,histname,bgname):
  """Draws a 4x2 grid with jT and Bg on top row for 4 jet pT bins and ratio between these on bottom row
  
  Args:
    dataset: Dataset to be used
    histname: Name of histogram to be used
    bgname: Name of background histogram to be used
  """
  jtHistos = dataset.getHist(histname,jetpt = True) 
  bgHistos = dataset.getHist(bgname,jetpt = True,isBg = True)
  fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d,shareY = False)
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
  
def JtWithBackgroundRatioAll(dataset,histname,bgname):
  """Draws a 4x4 grid with jT and Bg on top and third for 8 jet pT bins and ratio between these on second and bottom row
  
  Args:
    dataset: Dataset to be used
    histname: Name of histogram to be used
    bgname: Name of background histogram to be used
  """
  jtHistos = dataset.getHist(histname,jetpt = True) 
  bgHistos = dataset.getHist(bgname,jetpt = True,isBg = True)
  fig, axs = defs.makegrid(4,4,xlog=True,ylog=True,d=d,shareY = False)
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
  Mixed_FullJets_R04 = datasetMixed("FullR04",NFIN=0,range=5,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=5)
 

  #measJt,jetpt = MB_FullJets_R04.getHist('JetConeJtWeightBin',jetpt = True)
  
  compareSetsWithRatio((MB_FullJets_R05,MB_FullJets_R04),'JetConeJtWeightBin')
  #compareSetsWithRatio((triggered_FullJets_R05,triggered_FullJets_R04),'JetConeJtWeightBin')
  #Mixed_FullJets_R04.printStats()
  #MB_FullJets_R04.printStats()
  #compareSets((Mixed_FullJets_R04,MB_FullJets_R04),'JetConeJtWeightBin')

  #JtWithBackgroundRatioAll(Mixed_FullJets_R04, 'JetConeJtWeightBin', 'BgJtWeightBin')


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
    
        
