import logging
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
from ROOT import TF1,TMath,Double
from matplotlib import container
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
import numpy as np

d = dict(
    jetalg = r'Anti-$k_\mathrm{T}$, R=0.4',
    jettype = 'Full Jets',
    system = r'ALICE p-Pb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' ,
    trigger = 'kINT7/kEMCEJE',
    cut = r'$\left| \eta_\mathrm{jet} \right| < 0.25$'
    
)

labelsize= 15
styles = [23,24,25,27,30]
colors = [1,2,3,4,7]
Rebin = 2
boxwidth = 2.5
xlow = 0.1
xhigh = 10
mSize = 6

def main(): 
  
  JetPtBins = [5,10,20,30,40,60,80,100,150,500]
  jetPt = [(JetPtBins[i],JetPtBins[i+1]) for i in range(8)]
  print(jetPt)
  #JetPtCenter = [7.5,15,25,35,50,70,90,125,325]
  JetPtCenter = [6.5,12.45,23.34,33.83,46.75,67.73,88.01,116.11,194.61]
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
  iS = 0

  f = root_open("errors_test.root", 'read')

  gGausRMS = f.Get("gGausRMS{:02d}".format(iS))
  gGausRMSerr = f.Get("gGausRMS{:02d}_Systematics".format(iS))
  gGausYield = f.Get("gGausYield{:02d}".format(iS))
  gGausYielderr = f.Get("gGausYield{:02d}_Systematics".format(iS))
  gGammaRMS = f.Get("gGammaRMS{:02d}".format(iS))
  gGammaRMSerr = f.Get("gGammaRMS{:02d}_Systematics".format(iS))
  gGammaYield = f.Get("gGammaYield{:02d}".format(iS))
  gGammaYielderr = f.Get("gGammaYield{:02d}_Systematics".format(iS))
  
  start = 4
  iS = 0
  stats = [None if ij < start else f.Get("JetConeJtWeightBinNFin{:02d}JetPt{:02d}_Statistics".format(iS,ij)) for ij in range(8)]
  fits = [None if ij < start else f.Get("JetConeJtWeightBinNFin{:02d}JetPt{:02d}_FitFunction".format(iS,ij)) for ij in range(8)]
  print(fits)


  print(stats)

  n_figs = 2
  
  
  
#   if(n_figs == 2):
#     fig,axs = defs.makeRatio(xlog=True,ylog=True,d=d,shareY=False,figsize = (5,6),grid=False)
#   else:
#     fig, axs = defs.makegrid(n_figs/2,2,xlog=True,ylog=True,d=d,shareY=False,figsize= (10,7.5) if n_figs == 4 else (n_figs*15/8,7.5) )
#   axs = axs.reshape(n_figs)
#   if(n_figs == 2):
#     pT = jetPt[start]
#     print(pT)
#     axs[0].text(0.8,7,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'] + '\n' + r'${:02d}\:\mathrm{{GeV}}/c < p_{{\mathrm{{T,jet}}}} < {:02d}\:\mathrm{{GeV}}/c$'.format(pT[0],pT[1]),fontsize = 10)
#   else:
#     axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'],fontsize = 11)
#   
  ratios = []
  xs2 = []

  for jT,pT,ij,fit in zip(stats[start:],jetPt[start:],range(start,9),fits[start:]):
    color = colors[1]
    fig,axs = defs.makeRatio(xlog=True,ylog=True,d=d,shareY=False,figsize = (5,6),grid=False)
    axs = axs.reshape(n_figs)
    axs[0].text(0.8,7,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'] + '\n' + r'${:02d}\:\mathrm{{GeV}}/c < p_{{\mathrm{{T,jet}}}} < {:02d}\:\mathrm{{GeV}}/c$'.format(pT[0],pT[1]),fontsize = 10)
    ax = axs[0]
    xs = np.arange(0,xhigh,0.01).tolist()
    for ii in range(6):
      print(fit.GetParameter(ii))
    #B2 is Gauss normalization
    #B3 is Gamma normalization
    gauss = fit.Clone()
    gauss.SetParameter(3,0)
    gamma = fit.Clone()
    gamma.SetParameter(0,0)
    ys = [fit.Eval(x) for x in xs]
    ys2 = [gauss.Eval(x) for x in xs]
    ys3 = [gamma.Eval(x) for x in xs]
    ax.plot(xs,ys2,'b:',label="Narrow")
    ax.plot(xs,ys3,'r--',label="Wide")
    ax.plot(xs,ys,'k',label="Total")

    jT.SetMarkerColor(color)
    jT.SetMarkerStyle(24)
    jT.SetLineColor(color)
    jT.SetMarkerSize(mSize)
    plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label="ALICE",fmt='o',fillstyle='none') #Plot jT histogram, 
    line = plot.get_children()[0]
    line.set_markersize(mSize)
    if(True):
      line.set_markerfacecolor('none')
      #line.set_markeredgecolor(color)
      line.set_color(color)
#     line.set_markerfacecolor('none')
#     line.set_markeredgecolor(colors[c])
#     line.set_markerfacecolor('none')
#     line.set_markeredgecolor('none')
#     line.set_drawstyle('default')
#     line.set_linestyle('dashed')  
#     line.set_color(colors[c])
    if(n_figs > 2):
      ax.text(0.5,1e2,r'${:02d}\:\mathrm{{GeV}} < p_{{\mathrm{{T,jet}}}} < {:02d}\:\mathrm{{GeV}}$'.format(pT[0],pT[1])) 
    ax.set_xlim([xlow,xhigh]) #Set x-axis limits
    ax.set_ylim([1e-6,2e3]) #Set y-axis limits
    #ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
    x_ = Double()
    y1 = Double()
    xe = Double()
    ye = Double()
    NC = jT.GetN()
    y2 = []
    y2e = []
    xs=[]
    ex = []
    for ii in range(NC):
      jT.GetPoint(ii,x_,y1)
      xe = jT.GetErrorX(ii)
      ye = jT.GetErrorY(ii)
      x1 = x_*1.0
      xs.append(x1)
      ex.append(xe)
      if(y1 > 0):
        y2.append(y1/fit.Eval(x1))
        y2e.append(ye/fit.Eval(x1))
      else:
        y2.append(0)
        y2e.append(0)
    #ratio = jT.Clone()
    ratio= Graph(NC)

    print(xs[::5])
    print(y2[::5])
    print(ex[::5])
    print(y2e[::5])
    for x0,y0,x0e,y0e,i in zip(xs,y2,ex,y2e,range(NC)):
      ratio.SetPoint(i,x0,y0)
      ratio.SetPointError(i,x0e,x0e,y0e,y0e)
    ratios.append(ratio)
    print(ratio)
    ax = axs[1]
  #for ratio,pT,ax,color in zip(ratios,jetPt[start:],axs[n_figs/2:n_figs+1],colors[1:]):
    
    print("Debug")
    ratio.SetMarkerColor(color)
    ratio.SetLineColor(color)
    ratio.SetMarkerStyle(24)
    ax.plot([0,20],[1,1],'k--')
    #ratio.SetMarkerSize(mSize)
    plot = rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label=r"Ratio",fmt='o',fillstyle='none') #Plot jT histogram, 
    line = plot.get_children()[0]
    line.set_markersize(mSize)
    if(True):
      line.set_markerfacecolor('none')
      #line.set_markeredgecolor(color)
      line.set_color(color)
    #ax.plot(xs2,ratio)
    ax.set_yscale('linear')
    ax.set_xlim([xlow,xhigh]) #Set x-axis limits
    ax.set_ylim([0,2.75]) #Set y-axis limits  
    ax.set_ylabel('Signal/Fit',fontsize=14) #Add y-axis labels to left- and righmost subfigures

    
    handles, labels = axs[0].get_legend_handles_labels()
    handles = [handles[3],handles[2],handles[0],handles[1]]
    labels = [labels[3],labels[2],labels[0],labels[1]]
  
    handles = [container.ErrorbarContainer(h,has_xerr=False,has_yerr=True) if isinstance(h, container.ErrorbarContainer) else h for h in handles]
    axs[0].legend(handles,labels,loc = 'lower left',numpoints=1)
  
    axs[0].text(0.11,3e2,"ALICE",weight='bold')
  
    fig.align_labels()
    print("Save PythonFigures/JtSignalFinalFitJetPt{}.pdf".format(ij))
    plt.savefig("PythonFigures/JtSignalFinalFitJetPt{}.pdf".format(ij),format='pdf') #Save figure
    plt.show() #Draw figure on screen


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
  if('pythiaW' in kwargs):
    pythiaW = kwargs.get('pythiaW')
    for x,c,l in zip(pythiaW,("black","red","blue","green"),('-', '--', '-.', ':')):
      hPythiaW = x[0]
      hPythiaW.SetMarkerColor(c)
      hPythiaW.SetLineColor(c)
      hPythiaW.SetMarkerStyle(25)
      print(x[1])
      plot = rplt.errorbar(hPythiaW,xerr=False,yerr=False,emptybins=False,fillstyle='none',axes=ax)
      line1, = plt.plot([1,2,3], label=x[1], linestyle=l,color=c)
      line = plot.get_children()[0]
      line.set_markerfacecolor('none')
      line.set_markeredgecolor('none')
      line.set_drawstyle('default')
      line.set_linestyle(l)  
      line.set_color(c)
    
  if('pythiaN' in kwargs):
    pythiaN = kwargs.get('pythiaN')
    for x,c,l in zip(pythiaN,("black","red","blue","green"),('-', '--', '-.', ':')):
      hPythiaN = x[0]
      hPythiaN.SetMarkerColor(c)
      hPythiaN.SetLineColor(c)
      hPythiaN.SetMarkerStyle(24)
      plot = rplt.errorbar(hPythiaN,xerr=False,yerr=False,emptybins=False,fillstyle='none',axes=ax,label="_{}".format(x[1]))
      line = plot.get_children()[0]
      line.set_markerfacecolor('none')
      line.set_markeredgecolor('none')  
      line.set_drawstyle('default')
      line.set_linestyle(l)  
      line.set_color(c)  
  
  
  if('pythiaN' in kwargs):
    h2.SetMarkerColor('red')
    h2.SetLineColor('red')
    h2.SetMarkerStyle(25)
    h.SetMarkerStyle(24)
    plot1 = rplt.errorbar(h2,xerr=False,emptybins=False,axes=ax,label='ALICE wide')
    plot2 = rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,label = 'ALICE Narrow')
    line = plot1.get_children()[0]
    #line.set_markerfacecolor('none')
    line.set_markeredgecolor('red')     
    line = plot2.get_children()[0]
    #line.set_markerfacecolor('none')
    line.set_markeredgecolor('black')    
  elif('jussiW' in kwargs):
    h2.SetMarkerColor('red')
    h2.SetLineColor('red')
    h2.SetMarkerStyle(25)
    h.SetMarkerStyle(24)
    plot1 = rplt.errorbar(h2,xerr=False,emptybins=False,axes=ax,label=r'Jet $j_\mathrm{T}$ Wide')
    plot2 = rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,label=r'Jet $j_\mathrm{T}$ Narrow') 
    line = plot1.get_children()[0]
    #line.set_markerfacecolor('none')
    line.set_markeredgecolor('red')     
    line = plot2.get_children()[0]
    #line.set_markerfacecolor('none')
    line.set_markeredgecolor('black')       
  else:
    h2.SetMarkerColor('red')
    h2.SetLineColor('red')
    h2.SetMarkerStyle(25)
    h.SetMarkerStyle(24)
    plot1 = rplt.errorbar(h2,xerr=False,emptybins=False,axes=ax,label='Wide')
    plot2 = rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,label='Narrow') 
    line = plot1.get_children()[0]
    #line.set_markerfacecolor('none')
    line.set_markeredgecolor('red')     
    line = plot2.get_children()[0]
    #line.set_markerfacecolor('none')
    line.set_markeredgecolor('black')
      


  if('jussiW' in kwargs):
    hJussiW = kwargs.get('jussiW')
    hJussiW.SetMarkerColor('red')
    hJussiW.SetLineColor('red')
    hJussiW.SetMarkerStyle(25)
    plot = rplt.errorbar(hJussiW,xerr=True,emptybins=False,fillstyle='none',axes=ax,label=r'Dihadron $j_\mathrm{T}$ Wide')
    line = plot.get_children()[0]
    line.set_markerfacecolor('none')
    line.set_markeredgecolor('red')
  if('jussiN' in kwargs):
    hJussiN = kwargs.get('jussiN')
    hJussiN.SetMarkerColor('black')
    hJussiN.SetLineColor('black')
    hJussiN.SetMarkerStyle(24)
    plot = rplt.errorbar(hJussiN,xerr=True,emptybins=False,fillstyle='none',axes=ax,label=r'Dihadron $j_\mathrm{T}$ Narrow')
    line = plot.get_children()[0]
    line.set_markerfacecolor('none')
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

  pc = PatchCollection(errorboxes, facecolor='0.65', alpha=0.6,edgecolor='None')
  ax.add_collection(pc)
  pc2 = PatchCollection(errorboxes2, facecolor='0.65', alpha=0.6,edgecolor='None')
  ax.add_collection(pc2)

  print("({} - {})/9.0 + {} is {}".format(yhigh,ylow,ylow,(yhigh-ylow)/9.0+ylow))
  if('titleLoc' in kwargs):
    x_tit,y_tit = kwargs.get('titleLoc')
    ax.text(x_tit,y_tit,title + '\n' + d['system'] +'\n'+ d['jettype'] + '\n' + d['jetalg'] + '\n' + d['cut'],
      fontsize = 10)
    
  else:
    if('jussiW' not in kwargs):  
      ax.text((xhigh-xlow)*0.1+xlow,4*(yhigh-ylow)/5.0+ylow,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'],
            fontsize = 10)
    else:
      ax.text((xhigh-xlow)*0.1+xlow,4*(yhigh-ylow)/5.0+ylow,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'],
            fontsize = 10)
  
  handles, labels = ax.get_legend_handles_labels()
  if(len(labels) > 4):
    ax.legend(np.delete(handles,4), np.delete(labels,4), loc = 'upper left')
    ax.text((xhigh-xlow)*0.09+xlow,0.57*(yhigh-ylow)+ylow,"ALICE",weight='bold')
  else:
    ax.legend(loc='upper right')
    if(len(labels) > 2):
      ax.text((xhigh-xlow)*0.770+xlow,0.66*(yhigh-ylow)+ylow,"ALICE",weight='bold')
    else:
      ax.text((xhigh-xlow)*0.770+xlow,0.80*(yhigh-ylow)+ylow,"ALICE",weight='bold')
  #ax.legend(loc = 'upper left') 
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  ax.tick_params(which='both',direction='in') #Move ticks from outside to inside

  plt.savefig("{}.pdf".format(file),format='pdf') #Save figure

  plt.show() #Draw figure on screen

   

   
if __name__ == "__main__": main()

    