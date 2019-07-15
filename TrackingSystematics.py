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
import ROOT
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
from dataset import *
labelsize= 15

def main(): 
  Rebin = 4
  colors = [7,1,2,4,6]
  Pythia = dataset("EFf 100",NFIN=1,range=(2,8),filename="Pythia/pythia8226_pp5TeV_Monash2013_Rbinned.root",directory='/JCDijetBaseTask/jcdijet',color=2,style=24,rebin=Rebin)
  Pythia2 = dataset("Eff 97",NFIN=1,range=(2,8),filename="Pythia/pythia8226_pp5TeV_Monash2013_eff97_14k.root",directory='/JCDijetBaseTask/jcdijet',color=3,style=24,rebin=Rebin)
  inclusive,jetPt = Pythia.getHist('JetConeJtWeightBin',jetpt = True)  
  datasets = [Pythia]
  datasets.append(Pythia2)
  incs = [inclusive]
  outputfile = "TrackingSyst.root"
  for data in datasets[1:]:
    incs.append(data.getHist('JetConeJtWeightBin',jetpt = False))
    
  names = [data.name() for data in datasets]
  signals  = [data.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = False) for data in datasets]
  graphs = [data.getGraphs() for data in datasets]
  gausRMS = [x[0] for x in graphs]
  gammaRMS = [x[1] for x  in graphs]
  gausYield = [x[2] for x in graphs]
  gammaYield = [x[3] for x in graphs]
  systGausRMS = defs.makeSystError(gausRMS[0],gausRMS[1])
  systGausRMS.SetName("GausRMS_trackingSystematics")
  systGammaRMS = defs.makeSystError(gammaRMS[0],gammaRMS[1])
  systGammaRMS.SetName("GammaRMS_trackingSystematics")
  systGausRMS_abs = defs.makeSystError(gausRMS[0],gausRMS[1],abs=True)
  systGausRMS_abs.SetName("GausRMS_trackingSystematics_absolute")
  systGammaRMS_abs = defs.makeSystError(gammaRMS[0],gammaRMS[1],abs=True)
  systGammaRMS_abs.SetName("GammaRMS_trackingSystematics_absolute")
  systGausRMS_rel = defs.makeSystError(gausRMS[0],gausRMS[1],rel=True)
  systGausRMS_rel.SetName("GausRMS_trackingSystematics_relative")
  systGammaRMS_rel = defs.makeSystError(gammaRMS[0],gammaRMS[1],rel=True)
  systGammaRMS_rel.SetName("GammaRMS_trackingSystematics_relative")
  defs.drawErrors2(systGausRMS_rel,40,150,0,-0.1,0.1,0,"Narrow RMS","Tracking","_data","SystematicErrors/SystematicErrorsGausRMS_Tracking.pdf",0,0,0.04);
  defs.drawErrors2(systGammaRMS_rel,40,150,0,-0.15,0.15,0,"Wide RMS","Tracking","_data","SystematicErrors/SystematicErrorsGammaRMS_Tracking.pdf",0,0,0.05);

  with root_open(outputfile,'recreate') as f:
    f.cd()
    systGausRMS.Write()
    systGammaRMS.Write()
    systGausRMS_abs.Write()
    systGammaRMS_abs.Write()
    systGausRMS_rel.Write()
    systGammaRMS_rel.Write()
  
  
  start = 2
  end = 7
  n_figs = 8
  n_rows = 2
  fig, axs = defs.makegrid(4,n_figs//4,xlog=True,ylog=True,d=d,shareY=False,figsize=(15,7.5))
  axs = axs.reshape(n_figs)
  axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone' + '\n Inclusive jT',fontsize = 7)
  for inc,name,color,j in zip(incs,names,colors,range(10)):
    print("Plot {}".format(name))
    for jT,pT,ax,i in zip(inc[start:],jetPt[start:],axs[0:4],range(0,9)):
      jT.SetMarkerColor(color)
      jT.SetLineColor(color)
      plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=name,fmt='o',fillstyle='none',ecolor='blue') #Plot jT histogram, 
      ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
  
      ax.set_xlim([0.1,22]) #Set x-axis limits
      ax.set_ylim([5e-4,2e3]) #Set y-axis limits
      ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
    ratios = []
    if j > 0:
      for jT,div in zip(inc,incs[0]):
        h = jT.Clone()
        h.Divide(div)
        ratios.append(h)
      axs[4].set_ylabel('Ratio',fontsize=9) #Add y-axis labels to left- and righmost subfigures
      axs[-1].set_ylabel('Ratio',fontsize=9)       
      for ratio,pT,ax in zip(ratios[start:],jetPt[start:],axs[4:9]):
        rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram,        
        #if(i == 0):
        ax.set_yscale('linear')
        ax.set_xlim([0.1,20]) #Set x-axis limits
        ax.set_ylim([0.8,1.2]) #Set y-axis limits     
  
    
  axs[0].legend(loc = 'lower left')
  plt.savefig("PythonFigures/TrackingR04JetConeJtInclusivePtFrom{}To{}.pdf".format(start,end),format='pdf') #Save figure
  plt.show() #Draw figure on screen

  fig, axs = defs.makegrid(4,n_figs//4,xlog=True,ylog=True,d=d,shareY=False,figsize=(15,7.5))
  axs = axs.reshape(n_figs)
  axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone' + '\n Subtracted jT',fontsize = 7)
  for signal,name,color,j in zip(signals,names,colors,range(10)):
    print("Plot {}".format(name))
    for jT,pT,ax,i in zip(signal[start:],jetPt[start:],axs[0:4],range(0,9)):
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
      axs[4].set_ylabel('Ratio',fontsize=9) #Add y-axis labels to left- and righmost subfigures
      axs[-1].set_ylabel('Ratio',fontsize=9)       
      for ratio,pT,ax in zip(ratios[start:],jetPt[start:],axs[4:9]):
        rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram,        
        #if(i == 0):
        ax.set_yscale('linear')
        ax.set_xlim([0.1,20]) #Set x-axis limits
        ax.set_ylim([0.8,1.2]) #Set y-axis limits     
  
    
  axs[0].legend(loc = 'lower left')
  plt.savefig("PythonFigures/TrackingR04JetConeJtSignalPtFrom{}To{}.pdf".format(start,end),format='pdf') #Save figure
  plt.show() #Draw figure on screen

  drawWithErrors2Combined(gausRMS,gammaRMS,15,500,1,0,1.65,0,r'jet $p_T$',r'$\sqrt{\left<j_T^2\right>}$','Pythia','PythonFigures/TrackingSystematicsRMS',separate=True)
  

def drawWithErrors2Combined(hists,hists2,xlow,xhigh,logx,ylow,yhigh,ylog,xTitle,yTitle,title,file,separate=False):
  colors = ['white','black','red','green','blue']

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
  names = ["100%","97%"]
  for h,h2,c,name in zip(hists,hists2,(1,2,3,4),names):
    ax = axs[0]
    h.SetMarkerColor(c)
    h.SetMarkerStyle(24)
    h2.SetMarkerColor(c)
    h2.SetMarkerStyle(25)
    ratios1.append(grrDivide(h, hists[1]))
    ratios2.append(grrDivide(h2, hists2[1]))
    if(separate):
      plot = rplt.errorbar(h,xerr=False,emptybins=False,label = name,axes=ax,fmt='+')
    else:
      if c == 1:
        plot = rplt.errorbar(h,xerr=False,emptybins=False,label = 'Narrow',axes=ax,fmt='+')
      else:
        plot = rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,fmt='+')
    #rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,fmt='+')
    line = plot.get_children()[0]
    #line.set_markerfacecolor('none')
    #line.set_markeredgecolor(colors[c])
#     line.set_markerfacecolor('none')
#     line.set_markeredgecolor(colors[c])
#     line.set_markerfacecolor('none')
#     line.set_markeredgecolor('none')
#     line.set_drawstyle('default')
#     line.set_linestyle('dashed')  
#     line.set_color(colors[c])
    if(separate):
      ax.set_xlim([xlow,xhigh])
      ax.set_ylim([ylow,yhigh])
      ax = axs[1]
      ax.grid(True)
      ax.set_xlim([xlow,xhigh])
      ax.set_ylim([ylow,yhigh])
    if(c > 1):
      if(separate):
        plot = rplt.errorbar(h2,xerr=False,yerr=True,emptybins=False,axes=ax,label=name,fmt='+')
      else:
        plot = rplt.errorbar(h2,xerr=False,yerr=True,emptybins=False,axes=ax,label='{} Wide'.format(name),fmt='+')
#       line = plot.get_children()[0]
#       line.set_markerfacecolor('none')
#       line.set_markeredgecolor('none')
#       line.set_drawstyle('default')
#       line.set_linestyle('solid')  
#       line.set_color(colors[c])  
    else: 
      rplt.errorbar(h2,xerr=False,emptybins=False,axes=ax,label='{}'.format(name),fmt='+')

  
  print("({} - {})/9.0 + {} is {}".format(yhigh,ylow,ylow,(yhigh-ylow)/9.0+ylow))
  if(separate):
    axs[1].text(65,0.2,title + '\n' + d['system'] +'\n'+  d['jettype'] + '\n' + r'Anti-$k_T$',fontsize = 10)
    axs[1].text(102,1.22,"Wide")
    axs[0].text(102,1.22,"Narrow")

  else:
    axs[0].text(100,(yhigh-ylow)/3.0+ylow,title + '\n' + d['system'] +'\n'+  d['jettype'] + '\n' + r'Anti-$k_T$',fontsize = 10)
   
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

  for ratio,ratio2,c in zip(ratios1,ratios2,(1,2,3,4)):
    #ratio.Shift(-2)
    #ratio2.Shift(2)
    ratio.SetMarkerColor(c)
    ratio2.SetMarkerColor(c)
    ratio.SetMarkerStyle(24)
    ratio2.SetMarkerStyle(25)
    plot = rplt.errorbar(ratio,xerr = False,yerr=True,emptybins=False,axes=ax)
#     line = plot.get_children()[0]
#     line.set_markerfacecolor('none')
#     line.set_markeredgecolor(colors[c])
#     line.set_markerfacecolor('none')
#     line.set_markeredgecolor('none')
#     line.set_drawstyle('default')
#     line.set_linestyle('dashed')  
#     line.set_color(colors[c])
    
    plot = rplt.errorbar(ratio2,xerr = False,yerr=True,emptybins=False,axes=ax2)
        
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

def grrDivide(gr1,gr2):
  x1 = ROOT.Double()
  x_ = ROOT.Double()
  x1e = ROOT.Double()
  x2 = ROOT.Double()
  x2e = ROOT.Double()
  y1 = ROOT.Double()
  y1e = ROOT.Double()
  y2 = ROOT.Double()
  y2e = ROOT.Double()
  test = []
  y = []
  ex = []
  ey = []

  NC =  gr1.GetN()
  x = [0 for i in range(NC)]
  test = [0 for i in range(NC)]

  for ii in range(NC):
    x1e = gr1.GetErrorX(ii)
    y1e = gr1.GetErrorY(ii)
    gr1.GetPoint(ii, x_, y1)
    x1 = x_ * 1.0
    gr2.GetPoint(ii, x2, y2)
    x2e = gr2.GetErrorX(ii)
    y2e = gr2.GetErrorY(ii)
    x[ii] = x1
    y.append(y1/y2 if y2 !=0 else 0)
    ex.append(x2e)
    ey.append(math.sqrt(pow(y1e/y2,2)+pow(y1*y2e/(y2*y2),2)) if y2 != 0 else 0)
    
  gr= Graph(NC)
  for x0,y0,x0e,y0e,i in zip(x,y,ex,ey,range(NC)):
    gr.SetPoint(i,x0,y0)
    gr.SetPointError(i,x0e,x0e,y0e,y0e)
  return gr
  
if __name__ == "__main__": main()
