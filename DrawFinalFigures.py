import logging
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
import rootpy
import defs
import re
import matplotlib
from rootpy.plotting import Canvas
from rootpy.plotting import Legend
from rootpy.plotting import Graph
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import math
from dataset import *

d = dict(
    jetalg = r'Anti-$k_\mathrm{T}$, R=0.4',
    jettype = 'Full Jets',
    system = r'ALICE pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' ,
    trigger = 'kINT7/kEMCEJE'
)

labelsize= 15
styles = [20,24,25,27,30]
colors = [1,2,3,4,7]
Rebin = 2
boxwidth = 2.5

def main(): 
  
  JetPtBins = [5,10,20,30,40,60,80,100,150,500]
  JetPtCenter = [7.5,15,25,35,50,70,90,125,325]
  LeadPtMode = [2.0,3.0,5.6,8.9,9.9,14.8,20.7,26.0,34.6]
  LeadPtMean = [2.35,3.72,6.66,9.59,13.58,17.98,23.27,27.55,31.68]
  LeadPtError = [0.93,1.69,3.18,4.43,6.74,8.34,9.68,10.27,10.55]
  JetPtError = [2.5,5,5,5,10,10,10,25,175]
  JetPtLeadPtG = Graph(len(JetPtCenter)) #Gives jet pT as a function of leading pT
  JetPtLeadPtGerr = Graph(len(JetPtCenter))
  LeadPtJetPtG = Graph(len(JetPtCenter))
  LeadPtJetPtGerr = Graph(len(JetPtCenter))
  for i,(x,y,e) in enumerate(zip(LeadPtMean,JetPtCenter,JetPtError)):
    JetPtLeadPtG.SetPoint(i,x,y)
    JetPtLeadPtGerr.SetPoint(i,x,e)
  for i,(x,y,e) in enumerate(zip(JetPtCenter,LeadPtMean,LeadPtError)):
    LeadPtJetPtG.SetPoint(i,x,y)
    LeadPtJetPtGerr.SetPoint(i,x,e)
  
  Njets = 8
  

  
  
  topcomment = "_systematics_Triggered"
  
  finderName = ["Full jets, Anti-k_\mathrm{T} R = 0.4","Charged Jets, Anti-k_\mathrm{T} R = 0.4"]
  finderType = ["Full","Charged"]
  setTitle = ["0","1","2","3","4","5"]
  finderR = {4,4}
  
  iC = 0
  logx = 1
  doWeight = 1
  mSize = 0.3
  iS = 0

  f = root_open("errors.root", 'read')
  f2 = root_open("jtAnalysisFinalGraphs.root",'read')

  gGausRMS = f.Get("gGausRMS{:02d}".format(iS))
  gGausRMSerr = f.Get("gGausRMS{:02d}_Systematics".format(iS))
  gGausYield = f.Get("gGausYield{:02d}".format(iS))
  gGausYielderr = f.Get("gGausYield{:02d}_Systematics".format(iS))
  gGammaRMS = f.Get("gGammaRMS{:02d}".format(iS))
  gGammaRMSerr = f.Get("gGammaRMS{:02d}_Systematics".format(iS))
  gGammaYield = f.Get("gGammaYield{:02d}".format(iS))
  gGammaYielderr = f.Get("gGammaYield{:02d}_Systematics".format(iS))


  

  g = gGausRMSerr
  n = g.GetN()
  xs = g.GetX()
  ys = g.GetY()
  xerrs = g.GetEX()
  yerrs = g.GetEY()
  for i,x,y,xerr,yerr in zip(range(n),xs,ys,xerrs,yerrs):
    if(y > 0 and yerr/y > 0.13):
      yerr = 0.13*y
    if(y > 0 and yerr/y < 0.09):
      yerr = 0.09*y
    g.SetPointError(i,xerr,yerr)
    
  g = gGammaRMSerr
  n = g.GetN()
  xs = g.GetX()
  ys = g.GetY()
  xerrs = g.GetEX()
  yerrs = g.GetEY()
  for i,x,y,xerr,yerr in zip(range(n),xs,ys,xerrs,yerrs):
    if(y > 0 and yerr/y > 0.13):
      yerr = 0.13*y
    if(y > 0 and yerr/y < 0.10):
      yerr = 0.08*y
    g.SetPointError(i,xerr,yerr)
    
  #gGausRMS_trackingerr = f2.Get("....")
  #gGammaRMS_trackingerr = f2.Get("....")
  #gGausRMSerr = defs.combineErrors(gGausRMSerr,gGausRMS_trackingerr)
  #gGammaRMSerr = defs.combineErrors(gGammaRMSerr,gGammaRMS_trackingerr)
  #gGausRMSerr = defs.addConstantError(gGausRMSerr,0.20,rel=True)
  #gGammaRMSerr = defs.addConstantError(gGammaRMSerr,0.09,rel=True)
  iSPythia = 3
  gGausRMSPythia = f.Get("gGausRMS{:02d}".format(iSPythia))
  gGausYieldPythia = f.Get("gGausYield{:02d}".format(iSPythia))
  gGammaRMSPythia = f.Get("gGammaRMS{:02d}".format(iSPythia))
  gGammaYieldPythia = f.Get("gGammaYield{:02d}".format(iSPythia))

  gGausRMSJussi = f2.Get("systematicError_RMSNarrow_p-Pb_5.02TeV_0.2<xlong<0.4")
  gGammaRMSJussi = f2.Get("systematicError_RMSWide_p-Pb_5.02TeV_0.2<xlong<0.4")

  doWhat = False
  n = gGausRMS.GetN()
  xs = gGausRMSerr.GetX()

  if(doWhat):
    for g in (gGausRMS,gGausRMSerr,gGausYield,gGausYielderr,gGammaRMS,gGammaRMSerr,gGammaYield,gGammaYielderr,gGausRMSPythia,gGausYieldPythia,gGammaRMSPythia,gGammaYieldPythia):
      n = g.GetN()
      xs = g.GetX()
      ys = g.GetY()
      xerrs = g.GetEX()
      yerrs = g.GetEY()
      for i,x,y,xerr,yerr in zip(range(n),xs,ys,xerrs,yerrs):
        g.SetPoint(i,LeadPtMean[i],y)
        g.SetPointError(i,LeadPtError[i],yerr)
  else:      
    for g in (gGausRMSJussi,gGammaRMSJussi):
      n = g.GetN()
      xs = g.GetX()
      ys = g.GetY()
      xerrs = g.GetEX()
      yerrs = g.GetEY()
      for i,x,y,xerr,yerr in zip(range(n),xs,ys,xerrs,yerrs):
        print("Trigger pT: {}".format(x))
        x1 = JetPtLeadPtG.Eval(x)
        x1e = JetPtLeadPtGerr.Eval(x)
        print("Estimated jet pT: {}".format(x1))
        g.SetPoint(i,x1,y)
        g.SetPointError(i,x1e,yerr)


  Pythia = dataset("Pythia8 4C",NFIN=0,range=(1,8),filename="Pythia/Grid_Monash.root",directory='/JCDijetBaseTask/jcdijet',color=colors[1],style=24,rebin=Rebin)
  Pythia2 = dataset("Pythia8 Monash",NFIN=0,range=(1,8),filename="Pythia/Grid_Tune4c.root",directory='/JCDijetBaseTask/jcdijet',color=colors[2],style=24,rebin=Rebin)
  Herwig = dataset("Herwig 7.0",NFIN=0,range=(1,8),filename="Herwig/Herwig-LHCtest.root",directory='/JJetJt',color=colors[3],style=24,rebin=Rebin)
  #Pythia_ALICE = dataset("ALICE Pythia6 Perugia2011",NFIN=0,range=(1,8),filename="CF_pPb_MC_legotrain/legotrain_610_20181010-1926_LHCb4_fix_CF_pPb_MC_ptHardMerged.root",directory='AliJJetJtTask/AliJJetJtHistManager',color=colors[4],style=24,rebin=Rebin)
  #datasets = [Pythia]
  #inclusive,jetPt = Pythia.getHist('JetConeJtWeightBin',jetpt = True)   
  #datasets.append(Mixed_FullJets_R04)
  inclusive,jetPt = Pythia.getHist('JetConeJtWeightBin',jetpt = True)
  datasets = [Pythia]
  incs = [inclusive]
  datasets.append(Pythia2)
  datasets.append(Herwig)
  for data in datasets[1:]:
    incs.append(data.getHist('JetConeJtWeightBin',jetpt = False))
    
  names = [data.name() for data in datasets]
  signals  = [data.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = False) for data in datasets]
  graphs = [data.getGraphs() for data in datasets]
  gausRMS = [x[0] for x in graphs]
  gammaRMS = [x[1] for x  in graphs]
  gausYield = [x[2] for x in graphs]
  gammaYield = [x[3] for x in graphs]  


    
  drawWithErrors2(gGausRMS,gGausRMSerr,40,150,0,0,0.6,0,r'$p_\mathrm{T,jet}$',r'$\sqrt{\left<j_\mathrm{T}^{2}\right>}$',"Narrow Component","_triggered","PythonFigures/gausRMSWithSystematics",1,-1);
  #drawWithErrors2(gGausYield,gGausYielderr,40,150,0,0,9,0,r'$p_\mathrm{T,jet}$',"Yield","Narrow Component","_triggered","PythonFigures/gausYieldWithSystematics",1,-1);
  drawWithErrors2(gGammaRMS,gGammaRMSerr,40,150,0,0,1.5,0,r'$p_\mathrm{T,jet}$',r'$\sqrt{\left<j_\mathrm{T}^{2}\right>}$',"Wide Component","_triggered","PythonFigures/gammaRMSWithSystematics",1,-1);
  #drawWithErrors2(gGammaYield,gGammaYielderr,40,150,0,0,10,0,r'$p_\mathrm{T,jet}$',"Yield","Wide Component","_triggered","PythonFigures/gammaYieldWithSystematics",1,-1);


  pythiaN = [(gGausRMSPythia,"Pythia6"),(gausRMS[0],"Pythia8 4C"),(gausRMS[1],"Pythia8 Monash"),(gausRMS[2],"Herwig 7.0")]
  pythiaW = [(gGammaRMSPythia,"Pythia6"),(gammaRMS[0],"Pythia8 4C"),(gammaRMS[1],"Pythia8 Monash"),(gammaRMS[2],"Herwig 7.0")]
  
  drawWithErrors2Combined(gGausRMS,gGausRMSerr,gGammaRMS,gGammaRMSerr,40,150,0,0,1.5,0,r'$p_\mathrm{T,jet}$',r'$\sqrt{\left<j_\mathrm{T}^{2}\right>}$',"RMS","_triggered","PythonFigures/RMSWithSystematics",1,-1,wideE=0.041,narrowE=0.051)
  drawWithErrors2Combined(gGausRMS,gGausRMSerr,gGammaRMS,gGammaRMSerr,40,150,0,0,1.9,0,r'$p_\mathrm{T,jet}$',r'$\sqrt{\left<j_\mathrm{T}^{2}\right>}$',"","_triggered","PythonFigures/RMSWithSystematics_Pythia",1,-1,pythiaN=pythiaN,pythiaW=pythiaW,wideE=0.041,narrowE=0.051,titleLoc=(100,1.5))
  if(doWhat):
    drawWithErrors2Combined(gGausRMS,gGausRMSerr,gGammaRMS,gGammaRMSerr,0,50,0,0,2.5,0,r'$p_\mathrm{T,trigger}$',r'$\sqrt{\left<j_\mathrm{T}^{2}\right>}$',"","_triggered","PythonFigures/RMSWithSystematics_DihadronTriggerPt",1,-1,jussiN=gGausRMSJussi,jussiW=gGammaRMSJussi,wideE=0.041,narrowE=0.051)
  else:
    drawWithErrors2Combined(gGausRMS,gGausRMSerr,gGammaRMS,gGammaRMSerr,10,150,0,0,2.5,0,r'$p_\mathrm{T,jet}$',r'$\sqrt{\left<j_\mathrm{T}^{2}\right>}$',"","_triggered","PythonFigures/RMSWithSystematics_DihadronJetPt",1,-1,jussiN=gGausRMSJussi,jussiW=gGammaRMSJussi,wideE=0.041,narrowE=0.051)

  #compareToJussi(gGausRMS,gGausRMSerr,gGausRMSJussi,5,150,0,0,0.7,0,r'$p_\mathrm{T,jet}$',r'$\sqrt{\left<j_\mathrm{T}^{2}\right>}$',"RMS","_triggered","PythonFigures/RMSWithSystematics_Pythia",1,-1)
  #compareToJussi(gGammaRMS,gGammaRMSerr,gGammaRMSJussi,5,150,0,0,1.5,0,r'$p_\mathrm{T,jet}$',r'$\sqrt{\left<j_\mathrm{T}^{2}\right>}$',"RMS","_triggered","PythonFigures/RMSWithSystematics_Pythia",1,-1)

  #errGraph = [f.Get('JetConeJtWeightBinNFin{:02d}JetPt{:02d}_Systematics'.format(iS,i)) for i in range(2,Njets)]  #Get jT histograms from file an array
  #hJtSignalGraph = [f.Get('JetConeJtWeightBinNFin{:02d}JetPt{:02d}_Statistics'.format(iS,i)) for i in range(2,Njets)]  #Get jT histograms from file an array

  for i in range(2,Njets):
    ylow = 5e-6
    yhigh = 5e2
    if(logx):
      xlow = 0.1
      xhigh = 10
    else:
      xlow  = 0
      xhigh = 5   
  #drawWithErrors2(hJtSignalGraph, errGraph, xlow, xhigh, logx, ylow/1000, yhigh*100 , 1, "pPb #sqrt{s_{NN}} = 5.02 TeV",finderName[1],"_triggered","FinalFigures/JetConejtWithSystematics", 1);


def compareToJussi(h,h_sys,h_jussi,xlow,xhigh,logx,ylow,yhigh,ylog,xTitle,yTitle,title,comment,file,iS,ij):
  ax = plt.gca()
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  ax.set_xlabel(xTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  ax.set_ylabel(yTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  h_jussi.SetMarkerColor('red')
  h_jussi.SetLineColor('red')
  rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,label='jet jT',fmt='+') #Plot jT histogram, 
  rplt.errorbar(h_jussi,xerr=False,emptybins=False,axes=ax,label='diHadron jT')
  errorboxes = []
  n = h_sys.GetN()
  xs = h_sys.GetX()
  ys = h_sys.GetY()
  xerrs = h_sys.GetEX()
  yerrs = h_sys.GetEY()
  for x in (xs,ys,xerrs,yerrs):
    x.SetSize(n)
  print(xs)
  print(ys)
  print(xerrs)
  print(yerrs)
  for x, y, xe, ye in zip(xs, ys, xerrs, yerrs):
    rect = Rectangle((x - boxwidth, y - ye), boxwidth*2,ye*2)
    errorboxes.append(rect)

  pc = PatchCollection(errorboxes, facecolor='r', alpha=0.5,edgecolor='None')
  ax.add_collection(pc)
  print("({} - {})/9.0 + {} is {}".format(yhigh,ylow,ylow,(yhigh-ylow)/9.0+ylow))
  ax.text(100,(yhigh-ylow)/12.0+ylow,title + '\n' + d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'],
          fontsize = 10)
    
  ax.legend(loc = 'lower left')
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  plt.savefig("{}.pdf".format(file),format='pdf') #Save figure

  plt.show() #Draw figure on screen
    

def drawWithErrors2(h,h_sys,xlow,xhigh,logx,ylow,yhigh,ylog,xTitle,yTitle,title,comment,file,iS,ij):
  ax = plt.gca()
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  ax.set_xlabel(xTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  ax.set_ylabel(yTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,label='jT',fmt='+') #Plot jT histogram, 
  errorboxes = []
  n = h_sys.GetN()
  xs = h_sys.GetX()
  ys = h_sys.GetY()
  xerrs = h_sys.GetEX()
  yerrs = h_sys.GetEY()
  for x in (xs,ys,xerrs,yerrs):
    x.SetSize(n)
  print(xs)
  print(ys)
  print(xerrs)
  print(yerrs)
  for x, y, xe, ye in zip(xs, ys, xerrs, yerrs):
    rect = Rectangle((x - boxwidth, y - ye), boxwidth*2,ye*2)
    errorboxes.append(rect)

  pc = PatchCollection(errorboxes, facecolor='r', alpha=0.5,edgecolor='None')
  ax.add_collection(pc)
  print("({} - {})/9.0 + {} is {}".format(yhigh,ylow,ylow,(yhigh-ylow)/9.0+ylow))
  ax.text(100,(yhigh-ylow)/12.0+ylow,title + '\n' + d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'],
          fontsize = 10)
    
  #ax.legend(loc = 'lower left')
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  plt.savefig("{}.pdf".format(file),format='pdf') #Save figure

  plt.show() #Draw figure on screen


def drawWithErrors2Combined(h,h_sys,h2,h2_sys,xlow,xhigh,logx,ylow,yhigh,ylog,xTitle,yTitle,title,comment,file,iS,ij,**kwargs):
  print("DrawWithErrors2Combined")
  ax = plt.gca()
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  ax.set_xlabel(xTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  ax.set_ylabel(yTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  h.SetMarkerColor('black')
  h.SetLineColor('black')
  if('pythiaN' in kwargs):
    h2.SetMarkerColor('red')
    h2.SetLineColor('red')
    h2.SetMarkerStyle(25)
    plot1 = rplt.errorbar(h2,xerr=False,emptybins=False,axes=ax,label='ALICE wide')
    plot2 = rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,label = 'ALICE Narrow')
    line = plot1.get_children()[0]
    line.set_markerfacecolor('none')
    line.set_markeredgecolor('red')     
    line = plot2.get_children()[0]
    line.set_markerfacecolor('none')
    line.set_markeredgecolor('black')    
  else:
    h2.SetMarkerColor('red')
    h2.SetLineColor('red')
    h2.SetMarkerStyle(25)
    plot1 = rplt.errorbar(h2,xerr=False,emptybins=False,axes=ax,label='Wide')
    plot2 = rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,label='Narrow') 
    line = plot1.get_children()[0]
    line.set_markerfacecolor('none')
    line.set_markeredgecolor('red')     
    line = plot2.get_children()[0]
    line.set_markerfacecolor('none')
    line.set_markeredgecolor('black')
      
  if('pythiaW' in kwargs):
    pythiaW = kwargs.get('pythiaW')
    for x,c in zip(pythiaW,("black","red","blue","green")):
      hPythiaW = x[0]
      hPythiaW.SetMarkerColor(c)
      hPythiaW.SetLineColor(c)
      hPythiaW.SetMarkerStyle(25)
      print(x[1])
      plot = rplt.errorbar(hPythiaW,xerr=False,yerr=False,emptybins=False,fillstyle='none',axes=ax, label=x[1])
      line = plot.get_children()[0]
      line.set_markerfacecolor('none')
      line.set_markeredgecolor('none')
      line.set_drawstyle('default')
      line.set_linestyle('solid')  
      line.set_color(c)
    
  if('pythiaN' in kwargs):
    pythiaN = kwargs.get('pythiaN')
    for x,c in zip(pythiaN,("black","red","blue","green")):
      hPythiaN = x[0]
      hPythiaN.SetMarkerColor(c)
      hPythiaN.SetLineColor(c)
      hPythiaN.SetMarkerStyle(24)
      plot = rplt.errorbar(hPythiaN,xerr=False,yerr=False,emptybins=False,fillstyle='none',axes=ax,label="_{}".format(x[1]))
      line = plot.get_children()[0]
      line.set_markerfacecolor('none')
      line.set_markeredgecolor('none')  
      line.set_drawstyle('default')
      line.set_linestyle('solid')  
      line.set_color(c)

  if('jussiW' in kwargs):
    hJussiW = kwargs.get('jussiW')
    hJussiW.SetMarkerColor('red')
    hJussiW.SetLineColor('red')
    hJussiW.SetMarkerStyle(25)
    plot = rplt.errorbar(hJussiW,xerr=True,emptybins=False,fillstyle='none',axes=ax,label='Dihadron Wide')
    line = plot.get_children()[0]
    line.set_markerfacecolor('red')
    line.set_markeredgecolor('red')
  if('jussiN' in kwargs):
    hJussiN = kwargs.get('jussiN')
    hJussiN.SetMarkerColor('black')
    hJussiN.SetLineColor('black')
    hJussiN.SetMarkerStyle(24)
    plot = rplt.errorbar(hJussiN,xerr=True,emptybins=False,fillstyle='none',axes=ax,label='Dihadron Narrow')
    line = plot.get_children()[0]
    line.set_markerfacecolor('black')
    line.set_markeredgecolor('black')
  print("Errorboxes")
  errorboxes = []
  n = h_sys.GetN()
  xs = h_sys.GetX()
  ys = h_sys.GetY()
  xerrs = h_sys.GetEX()
  yerrs = h_sys.GetEY()
  for x in (xs,ys,xerrs,yerrs):
    x.SetSize(n)
  yeC = kwargs.get('narrowE',0)
  for x, y, xe, ye in zip(xs, ys, xerrs, yerrs):
    if(yeC > 0):
      ye2 = yeC * y
      ye = math.sqrt(ye**2 + ye2**2)
    rect = Rectangle((x - boxwidth, y - ye), boxwidth*2,ye*2)
    errorboxes.append(rect)

  errorboxes2 = []
  n = h2_sys.GetN()
  xs = h2_sys.GetX()
  ys = h2_sys.GetY()
  xerrs = h2_sys.GetEX()
  yerrs = h2_sys.GetEY()
  for x in (xs,ys,xerrs,yerrs):
    x.SetSize(n)
    
  yeC = kwargs.get('wideE',0)
  for x, y, xe, ye in zip(xs, ys, xerrs, yerrs):
    if(yeC > 0):
      ye2 = yeC * y
      ye = math.sqrt(ye**2 + ye2**2)
    rect = Rectangle((x - boxwidth, y - ye), boxwidth*2,ye*2)
    errorboxes2.append(rect)

  pc = PatchCollection(errorboxes, facecolor='r', alpha=0.5,edgecolor='None')
  ax.add_collection(pc)
  pc2 = PatchCollection(errorboxes2, facecolor='r', alpha=0.5,edgecolor='None')
  ax.add_collection(pc2)

  print("({} - {})/9.0 + {} is {}".format(yhigh,ylow,ylow,(yhigh-ylow)/9.0+ylow))
  if('titleLoc' in kwargs):
    x_tit,y_tit = kwargs.get('titleLoc')
    ax.text(x_tit,y_tit,title + '\n' + d['system'] +'\n'+ d['jettype'] + '\n' + d['jetalg'],
      fontsize = 10)
    
  else:
    if('jussiW' not in kwargs):  
      ax.text(100,(yhigh-ylow)/3.0+ylow,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'],
            fontsize = 10)
    else:
      ax.text(100,2*(yhigh-ylow)/3.0+ylow,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'],
            fontsize = 10)  
  ax.legend(loc = 'upper left') 
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  plt.savefig("{}.pdf".format(file),format='pdf') #Save figure

  plt.show() #Draw figure on screen

   

   
if __name__ == "__main__": main()

    