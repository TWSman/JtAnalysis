import logging
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
import rootpy
import defs
import re
import matplotlib
from matplotlib import container

from rootpy.plotting import Canvas
from rootpy.plotting import Legend
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
from dataset import *
import ROOT
import math
labelsize = 15
mSize = 6
names = []
styles = [20,24,25,27,30]
colors = [1,2,3,4,7]
xlow = 0.1

d = dict(
    jetalg = r'Anti-$k_\mathrm{T}$, R=0.4',
    jettype = 'Full Jets',
    system = r'ALICE p-Pb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' ,
    trigger = 'kINT7/kEMCEJE',
    cut = r'$\left| \eta_\mathrm{jet} \right| < 0.25$'
    
)

def main(): 
  Rebin = 4
  Mixed_FullJets_R04_MB = datasetMixed("Minimum Bias",NFIN=0,range=(1,5),filename="CF_pPb_legotrain/legotrain_CF_pPb_2274_20181219/legotrain_CF_pPb_2274_20181219_LHC13cde.root",directory='AliJJetJtTask/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=colors[0],style=24,rebin=Rebin)
  Mixed_FullJets_R04_HM_01 = datasetMixed("V0A, 0 -   1%",NFIN=0,range=(1,5),filename="CF_pPb_legotrain/legotrain_CF_pPb_2274_20181219/legotrain_CF_pPb_2274_20181219_LHC13cde.root",directory='AliJJetJtTask_Central01/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE_Central01/AliJJetJtHistManager',color=colors[1],style=24,rebin=Rebin)
  Mixed_FullJets_R04_HM_10 = datasetMixed("V0A, 0 -  10%",NFIN=0, range=(1,5),filename="CF_pPb_legotrain/legotrain_CF_pPb_2305_20190109/legotrain_CF_pPb_2305_20190109_LHC13bcde.root",directory='AliJJetJtTask_Central10/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE_Central10/AliJJetJtHistManager',color=colors[2],style=24,rebin=Rebin)
  Mixed_FullJets_R04_HM_0_1 = datasetMixed("V0A, 0 - 0.1%",NFIN=0, range=(1,5),filename="CF_pPb_legotrain/legotrain_CF_pPb_2305_20190109/legotrain_CF_pPb_2305_20190109_LHC13bcde.root",directory='AliJJetJtTask_Central0_1/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE_Central0_1/AliJJetJtHistManager',color=colors[3],style=24,rebin=Rebin)

  inclusive,jetPt = Mixed_FullJets_R04_MB.getHist('JetConeJtWeightBin',jetpt = True)
  datasets = [Mixed_FullJets_R04_MB]
  incs = [inclusive]
  datasets.append(Mixed_FullJets_R04_HM_0_1)
  datasets.append(Mixed_FullJets_R04_HM_01)
  datasets.append(Mixed_FullJets_R04_HM_10)

  for data in datasets[1:]:
    incs.append(data.getHist('JetConeJtWeightBin',jetpt = False))
  names = [data.name() for data in datasets]
  signals  = [data.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = False) for data in datasets]
  background = [data.getHist('BgJtWeightBin',jetpt=False) for data in datasets]
  graphs = [data.getGraphs() for data in datasets]
  gausRMS = [x[0] for x in graphs]
  gammaRMS = [x[1] for x  in graphs]
  gausYield = [x[2] for x in graphs]
  gammaYield = [x[3] for x in graphs]
  
  #drawWithErrors2Combined(gausRMS,gammaRMS,names,15,500,1,0,1.65,0,r'jet $p_T$',r'$\sqrt{\left<j_\mathrm{T}^2\right>}$','Pythia','PythonFigures/RMScomparisonHighM',separate=True)
  
  start = 4
  end = 8
  n_figs = 2
  if(end > start + n_figs/2):
    end = start + n_figs/2
  n_rows = 2
  if (False):
    if(n_figs == 2):
      fig,axs = defs.makeRatio(xlog=True,ylog=True,d=d,shareY=False,figsize = (5,7.5))
    else:
      fig, axs = defs.makegrid(n_figs/2,2,xlog=True,ylog=True,d=d,shareY=False,figsize= (10,7.5) if n_figs == 4 else (n_figs*15/8,7.5) )
    axs = axs.reshape(n_figs)
    axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'],fontsize = 11)
    for inc,name,color,j in zip(incs,names,colors,range(10)):
      print("Plot {}".format(name))
      for jT,pT,ax,i in zip(inc[start:],jetPt[start:],axs[0:n_figs/2],range(0,9)):
        jT.SetMarkerColor(color)
        jT.SetMarkerStyle(styles[j])
        jT.SetLineColor(color)
        plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=name,fmt='o',fillstyle='none') #Plot jT histogram, 
        line = plot.get_children()[0]
        line.set_markersize(mSize)
        if(styles[j] > 23):
          line.set_markerfacecolor('none')
          #line.set_markeredgecolor(color)
          line.set_color(color)
        ax.text(0.5,1e2,r'${:02d}\:\mathrm{{GeV}} < p_\mathrm{{T,jet}} < {:02d}\:\mathrm{{GeV}}$'.format(pT[0],pT[1])) 
    
        ax.set_xlim([xlow,10]) #Set x-axis limits
        ax.set_ylim([5e-5,2e3]) #Set y-axis limits
        #ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
      ratios = []
      if j > 0:
        for jT,div in zip(inc,incs[0]):
          h = jT.Clone()
          h.Divide(div)
          ratios.append(h)
        axs[n_figs/2].set_ylabel('Ratio to {}'.format(names[0]),fontsize=14) #Add y-axis labels to left- and righmost subfigures
        if(n_figs > 4): axs[-1].set_ylabel('Ratio to {}'.format(names[0]),fontsize=18)       
        for ratio,pT,ax in zip(ratios[start:],jetPt[start:],axs[n_figs/2:n_figs+1]):
          plot = rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram, 
          ax.plot([0,10],[1,1],'k--')
          line = plot.get_children()[0]
          line.set_markersize(mSize)
          if(styles[j] > 23):
            line.set_markerfacecolor('none')
          line.set_color(color)            
          #if(i == 0):
          ax.set_yscale('linear')
          ax.set_xlim([xlow,10]) #Set x-axis limits
          ax.set_ylim([0,2.2]) #Set y-axis limits     
    
      
    axs[0].legend(loc = 'lower left')
    plt.savefig("PythonFigures/HighMJetConeJtInclusivePtFrom{}To{}.pdf".format(start,end),format='pdf') #Save figure
    plt.show() #Draw figure on screen
  
  if(n_figs == 2):
    fig,axs = defs.makeRatio(xlog=True,ylog=True,d=d,shareY=False,figsize = (5,7.5),grid=False)
  else:
    fig, axs = defs.makegrid(n_figs/2,2,xlog=True,ylog=True,d=d,shareY=False,figsize= (10,7.5) if n_figs == 4 else (n_figs*15/8,7.5) )
  axs = axs.reshape(n_figs)
  if(n_figs == 2):
    pT = jetPt[start]
    print(pT)
    axs[0].text(0.8,7,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'] + '\n' + r'${:02d}\:\mathrm{{GeV}} < p_{{\mathrm{{T,jet}}}} < {:02d}\:\mathrm{{GeV}}$'.format(pT[0],pT[1]),fontsize = 11)
  else:
    axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'],fontsize = 11)
  for signal,name,color,j in zip(signals,names,colors,range(10)):
    print("Plot {}".format(name))
    for jT,pT,ax,i in zip(signal[start:],jetPt[start:],axs[0:n_figs/2],range(0,9)):
      jT.SetMarkerColor(color)
      jT.SetMarkerStyle(styles[j])
      jT.SetLineColor(color)
      plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=name,fmt='o',fillstyle='none') #Plot jT histogram, 
      line = plot.get_children()[0]
      line.set_markersize(mSize)
      if(styles[j] > 23):
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
      ax.set_xlim([xlow,10]) #Set x-axis limits
      ax.set_ylim([1e-5,2e3]) #Set y-axis limits
      #ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
    ratios = []
    if j > 0:
      for jT,div in zip(signal,signals[0]):
        h = jT.Clone()
        h.Divide(div)
        ratios.append(h)
      axs[n_figs/2].set_ylabel('Ratio to {}'.format(names[0]),fontsize=12) #Add y-axis labels to left- and righmost subfigures
      if(n_figs > 4): axs[-1].set_ylabel('Ratio to {}'.format(names[0]),fontsize=12)       
      for ratio,pT,ax in zip(ratios[start:],jetPt[start:],axs[n_figs/2:n_figs+1]):
        plot = rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram,
        ax.plot([0,10],[1,1],'k--')
     
        line = plot.get_children()[0]
        line.set_markersize(mSize)
        if(styles[j] > 23):
          line.set_markerfacecolor('none')
        line.set_color(color)  
        #if(i == 0):
        ax.set_yscale('linear')
        ax.set_xlim([xlow,10]) #Set x-axis limits
        ax.set_ylim([0.5,1.5]) #Set y-axis limits     
  
  handles, labels = axs[0].get_legend_handles_labels()
  handles = [container.ErrorbarContainer(h,has_xerr=False,has_yerr=True) if isinstance(h, container.ErrorbarContainer) else h for h in handles]
  axs[0].legend(handles,labels,loc = 'lower left',numpoints=1,prop={'family': 'monospace'})  
  #axs[0].legend(loc = 'lower left')
  axs[0].text(0.11,3e2,"ALICE",weight='bold')
  fig.align_labels()
  plt.savefig("PythonFigures/HighMJetConeJtSignalPtFrom{}To{}.pdf".format(start,end),format='pdf') #Save figure
  plt.show() #Draw figure on screen


  return 
  if(n_figs == 2):
    fig,axs = defs.makeRatio(xlog=True,ylog=True,d=d,shareY=False,figsize = (5,7.5))
  else:
    fig, axs = defs.makegrid(n_figs/2,2,xlog=True,ylog=True,d=d,shareY=False,figsize= (10,7.5) if n_figs == 4 else (n_figs*15/8,7.5) )
  axs = axs.reshape(n_figs)
  axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n' + d['cut'],fontsize = 11)
  for signal,name,color,j in zip(background,names,colors,range(10)):
    print("Plot {}".format(name))
    for jT,pT,ax,i in zip(signal[start:],jetPt[start:],axs[0:n_figs/2],range(0,9)):
      jT.SetMarkerColor(color)
      jT.SetMarkerStyle(styles[j])
      jT.SetLineColor(color)
      plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=name,fmt='o',fillstyle='none') #Plot jT histogram, 
      line = plot.get_children()[0]
      line.set_markersize(mSize)
      if(styles[j] > 23):
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
      ax.text(0.5,1e2,r'${:02d}\:\mathrm{{GeV}} < p_{{T,\mathrm{{jet}}}} < {:02d}\:\mathrm{{GeV}}$'.format(pT[0],pT[1])) 
      ax.set_xlim([xlow,22]) #Set x-axis limits
      ax.set_ylim([5e-4,2e3]) #Set y-axis limits
      #ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
    ratios = []
    if j > 0:
      for jT,div in zip(signal,background[0]):
        h = jT.Clone()
        h.Divide(div)
        ratios.append(h)
      axs[n_figs/2].set_ylabel('Ratio to {}'.format(names[0]),fontsize=18) #Add y-axis labels to left- and righmost subfigures
      if(n_figs > 4): axs[-1].set_ylabel('Ratio to {}'.format(names[0]),fontsize=18)       
      for ratio,pT,ax in zip(ratios[start:],jetPt[start:],axs[n_figs/2:n_figs+1]):
        plot = rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram,        
        line = plot.get_children()[0]
        line.set_markersize(mSize)
        if(styles[j] > 23):
          line.set_markerfacecolor('none')
        line.set_color(color)  
        #if(i == 0):
        ax.set_yscale('linear')
        ax.set_xlim([xlow,9]) #Set x-axis limits
        ax.set_ylim([0,2.2]) #Set y-axis limits     
  
  handles, labels = axs[0].get_legend_handles_labels()
  handles = [container.ErrorbarContainer(h,has_xerr=False,has_yerr=True) if isinstance(h, container.ErrorbarContainer) else h for h in handles]
  axs[0].legend(handles,labels,loc = 'lower left',numpoints=1)
  #axs[0].legend(loc = 'lower left')
  plt.savefig("PythonFigures/HighMJetConeJtBackgroundPtFrom{}To{}.pdf".format(start,end),format='pdf') #Save figure
  plt.show() #Draw figure on screen



def drawWithErrors2Combined(hists,hists2,names,xlow,xhigh,logx,ylow,yhigh,ylog,xTitle,yTitle,title,file,separate=False):
  #colors = ['white','black','red','green','blue']

  if(separate):
    fig,axs = plt.subplots(2,2,sharex=True,figsize=(7,7))
    axs = axs.reshape(4)
    axs[0].grid(True)
    for ax in axs:
      ax.yaxis.set_ticks_position('both')
      ax.tick_params(which='both',direction='in')
    for ax in axs[1::2]:
      ax.yaxis.set_label_position('right')
      ax.yaxis.tick_right()
      ax.tick_params(which='both',direction='in')
  else:
    fig,axs = plt.subplots(2,1,sharex=True,figsize= (7,7))
  ax = axs[0]
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  ax.set_xlabel(xTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  ax.set_ylabel(yTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  if(separate):
    axs[1].set_ylabel(yTitle,fontsize=labelsize)
  ax.set_ylabel(yTitle,fontsize=labelsize) #Add x-axis labels for bottom row

  ratios1 = []
  ratios2 = []
  if logx:
    ax.set_xscale('log')
  print(names)
  #names = ["Data","Pythia8 Monash2013","Pythia8 4C","Herwig","Pythia6 Perugia2011"]
  for h,h2,c,name,i in zip(hists,hists2,colors,names,range(10)):
    
    ax = axs[0]
    h.SetMarkerColor(c)
    h.SetMarkerStyle(styles[i])
    h2.SetMarkerColor(c)
    h2.SetMarkerStyle(styles[i])
    ratios1.append(defs.grrDivide(h, hists[1]))
    ratios2.append(defs.grrDivide(h2, hists2[1]))
    if(separate):
      plot = rplt.errorbar(h,xerr=False,emptybins=False,label = name,axes=ax,fmt='+')
    else:
      if c == 1:
        plot = rplt.errorbar(h,xerr=False,emptybins=False,label = 'Narrow',axes=ax,fmt='+')
      else:
        plot = rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,fmt='+')
    #rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,fmt='+')
    line = plot.get_children()[0]
    line.set_markersize(mSize)
    if(styles[i] > 23):
      line.set_markerfacecolor('none')
    if(separate):
      ax.set_xlim([xlow,xhigh])
      ax.set_ylim([ylow,yhigh])
      ax = axs[1]
      ax.grid(True)
      ax.set_xlim([xlow,xhigh])
      ax.set_ylim([ylow,yhigh])
    if(i > 0):
      if(separate):
        plot = rplt.errorbar(h2,xerr=False,yerr=True,emptybins=False,axes=ax,label=name,fmt='+')
      else:
        plot = rplt.errorbar(h2,xerr=False,yerr=True,emptybins=False,axes=ax,label='{} Wide'.format(name),fmt='+')
    else: 
      plot = rplt.errorbar(h2,xerr=False,emptybins=False,axes=ax,label='{}'.format(name),fmt='+')
    line = plot.get_children()[0]
    line.set_markersize(mSize)
    if(styles[i] > 23):
      line.set_markerfacecolor('none')
  
  print("({} - {})/9.0 + {} is {}".format(yhigh,ylow,ylow,(yhigh-ylow)/9.0+ylow))
  if(separate):
    axs[1].text(65,0.2,title + '\n' + d['system'] +'\n'+  d['jettype'] + '\n' + r'Anti-$k_\mathrm{T}$',fontsize = 10)
    axs[1].text(102,1.22,"Wide")
    axs[0].text(102,1.22,"Narrow")

  else:
    axs[0].text(100,(yhigh-ylow)/3.0+ylow,title + '\n' + d['system'] +'\n'+  d['jettype'] + '\n' + r'Anti-$k_\mathrm{T}$',fontsize = 10)
   
  axs[0].legend(loc = 'upper left')
  if(separate):
    for ax in axs[0:2]:
      ax.set_xlim([xlow,xhigh])
      ax.set_ylim([ylow,yhigh])
      ax.grid(True)
  else:
    ax.set_xlim([xlow,xhigh])
    ax.set_ylim([ylow,yhigh])
  
  if(separate):
    ax = axs[2]
    ax2 = axs[3]
  else:
    ax = axs[1]
    ax2 = axs[1]
  ax.grid(True)
  ax2.grid(True)
  ax.set_xlabel(xTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  ax2.set_xlabel(xTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  ax.set_ylabel('Ratio',fontsize=labelsize) #Add x-axis labels for bottom row
  ax2.set_ylabel('Ratio',fontsize=labelsize) #Add x-axis labels for bottom row

  for ratio,ratio2,c,i in zip(ratios1[1:],ratios2[1:],colors,range(1,10)):
    #ratio.Shift(-2)
    #ratio2.Shift(2)
    ratio.SetMarkerColor(c)
    ratio2.SetMarkerColor(c)
    ratio.SetMarkerStyle(styles[i])
    ratio2.SetMarkerStyle(styles[i])
    plot = rplt.errorbar(ratio,xerr = False,yerr=True,emptybins=False,axes=ax)
    line = plot.get_children()[0]
    line.set_markersize(mSize)
    if(styles[i] > 23):
      line.set_markerfacecolor('none')

    plot = rplt.errorbar(ratio2,xerr = False,yerr=True,emptybins=False,axes=ax2)
    line = plot.get_children()[0]
    line.set_markersize(mSize)
    if(styles[i] > 23):
      line.set_markerfacecolor('none')
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([0.8,1.22])
  if(separate):
    ax2.set_xlim([xlow,xhigh])
    ax2.set_ylim([0.8,1.22]) 
    ax2.grid(True) 
  if logx:
    ax.set_xscale('log')
    if(separate):
      ax2.set_xscale('log')
  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0  
  plt.savefig("{}.pdf".format(file),format='pdf') #Save figure

  plt.show() #Draw figure on screen
  
if __name__ == "__main__": main()
