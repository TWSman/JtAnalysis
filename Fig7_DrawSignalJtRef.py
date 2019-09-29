import logging
from matplotlib.pyplot import box
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
import rootpy
import defs
import re
import matplotlib
from rootpy.plotting import Canvas,Graph,Legend
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
#from PythiaComparison import drawWithErrors2Combined
from PythiaComparison import grrDivide
from dataset import *
from ROOT import TF1,TMath,TGraphErrors
import ROOT
from array import array
from ctypes import c_int
from matplotlib import container

Njets = 8
B1low = 0.1
B1start = [0.15,0.15,0.17,0.18,0.15,0.15,0.21,0.15,0.15]
B1high = 0.5
B2low = 40
B2high = 150
B2start = [50,50,66.20,67.90,50,50,75.89,77.62,77.62]
B3cut = 6
B3low = 0
B3high1 = 10
B3high2 = 25
B3start = [7,7,5.06,4.90,10,10,7,12.93,12.93]
B4low = 2.0
B4high = 15
B4start = [3,3,9.91,8.80,3,3,5.62,4.18,4.18]
B5low = 0.90
B5high = 4.5
B5start = [1.5,1.5,4.5,4.5,1.5,1.5,2.88,1.62,1.62]
peakstart = 0.5
peaklow = 0.3
peakhigh = 0.8
Rs = (0.3,0.4,0.5,-0.4)
styles = [23,24,25,27,30]

labelsize= 15
mSize = 6


def main(): 
  print ('Number of arguments: ', len(sys.argv), 'arguments.')
  print ('Argument list:',str(sys.argv))
  filename = "rootFiles/legotrain_CF_pPb_2305_20190109_LHC13bcde_minimal.root"
  separate = 0
  start = 2
  end = 6
  n_figs = end-start

  Mixed_FullJets_R04 = datasetMixed("Full jets R=0.4",NFIN=0,range=(1,5),filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=2)
  compareHistsWithRatio(Mixed_FullJets_R04,['JtWeightBin','JtWeightLeadingRefBin','JtWeightLeadingRefBin'],['Jet axis ref.','leading ref. (xlong 0.0-0.2)','leading ref. (xlong 0.2-0.4)'],step=1,start=1,extras=['','Xlong00','Xlong01'])
  plt.savefig("PythonFigures/JetVsLeadingRefConst.pdf",format='pdf') #Save figure
  plt.show()
  sets = compareHistsWithRatio(Mixed_FullJets_R04,['JetConeJtWeightBin','JetConeJtWeightLeadingRefBin','JetConeJtWeightLeadingRefBin','JetConeJtWeightLeadingRefBin','JetConeJtWeightLeadingRefBin'],['Jet axis ref.','leading ref.(xlong 0.0-0.2)','leading ref.(xlong 0.2-0.4)','leading ref.(xlong 0.4-0.6)','leading ref. (xlong 0.6-1.0)'],step=1,extras=['','Xlong00','Xlong01','Xlong02','Xlong03'])

  plt.savefig("PythonFigures/JetVsLeadingRefJetCone.pdf",format='pdf') #Save figure
  plt.show()
  
  JtJet = sets[0][0]
  JtLeadingxlong00 = sets[1][0]
  JtLeadingxlong01 = sets[2][0]
  JtLeadingxlong02 = sets[3][0]
  JtLeadingxlong03 = sets[4][0]
  JtLeading = [h.clone() for h in sets[1][0]]
  for h,s,s2,s3 in zip(JtLeading,sets[2][0],sets[3][0],sets[4][0]):
    h.Add(s,1)
    h.Add(s2,1)
    h.Add(s3,1)

  jetPt = sets[0][1]
  jetPtCenter = array('d',[(a+b)/2.0 for a,b in jetPt])
  jetPtErrors = array('d',[(b-a)/2.0 for a,b in jetPt])
  
  FullJets_fit = []
  FullJets_parameters = []
  FullJets_gausRMS = []
  FullJets_gammaRMS = []
  FullJets_gausYield = []
  FullJets_gammaYield = []
  for jT,title in zip((JtJet,JtLeading,JtLeadingxlong00,JtLeadingxlong01,JtLeadingxlong02),("Jet ref.","Leading ref","Leading ref.(xlong 0.0-0.2)","Leading ref.(xlong 0.2-0.4)","Leading ref.(xlong 0.4-0.6)","Leading ref.(xlong 0.6-1.0)")):
    gausRMS = []
    gammaRMS = []
    gausRMSe = []
    gammaRMSe = []
    gausYield = []
    gammaYield = []
    gausYielde = []
    gammaYielde = []
    fits = []
    parameters = []
    for h,i in zip(jT,range(Njets)):
      fit,d = defs.fitJtHisto(h,'',1,i,8,title,draw=False)
      fits.append(fit)
      parameters.append(d)
      gausRMS.append(d['gausRMS'])
      gausRMSe.append(d['gausRMSe'])
      gammaRMS.append(d['gammaRMS'])
      gammaRMSe.append(d['gammaRMSe'])
      gausYield.append(d['gausYield'])
      gausYielde.append(d['gausYielde'])
      gammaYield.append(d['gammaYield'])
      gammaYielde.append(d['gammaYielde'])
      
    gausRMSg = Graph(len(gausRMS)-2)
    gammaRMSg = Graph(len(gammaRMS)-2)
    gausYieldg = Graph(len(gausYield) -2)
    gammaYieldg = Graph(len(gammaYield)-2)
    for h,he,g in zip((gausYield,gammaYield),(gausYielde,gammaYielde),(gausYieldg,gammaYieldg)):
      for x,xe,a,e,i in zip(jetPtCenter[2:],jetPtErrors[2:],h[2:],he[2:],range(len(gausRMS)-2)):
        g.SetPoint(i,x,a)
        g.SetPointError(i,xe,xe,e,e)
    
    for a,b,c,d,e,f,i in zip(gausRMS[2:],gammaRMS[2:],gausRMSe[2:],gammaRMSe[2:],jetPtCenter[2:],jetPtErrors[2:],range(len(gausRMS)-2)):
      gausRMSg.SetPoint(i,e,a)
      gausRMSg.SetPointError(i,f,f,c,c)
      gammaRMSg.SetPoint(i,e,b)
      gammaRMSg.SetPointError(i,f,f,d,d)
    
    
    FullJets_gausRMS.append(gausRMSg)
    FullJets_gammaRMS.append(gammaRMSg)
    FullJets_gausYield.append(gausYieldg)
    FullJets_gammaYield.append(gammaYieldg)
    FullJets_fit.append(fits)
    FullJets_parameters.append(parameters)
  
  print(gausRMS[2:])
  print(gammaRMS[2:])
  print(jetPtCenter[2:])
  
  drawWithErrors2Combined(FullJets_gausRMS,FullJets_gammaRMS,["Jet ref.","Leading ref.","Leading ref.(xlong 0.0-0.2)","Leading ref.(xlong 0.2-0.4)","Leading ref.(xlong 0.4-0.6)"],15,500,1,0,1.85,0,r'jet $p_T$',r'$\sqrt{\left<j_T^2\right>}$','Pythia','PythonFigures/JetVsLeadingRefJetConeFits')
  

  if(separate > 0):
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.text(0.2,0.0005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone',fontsize = 7)
    rplt.errorbar(signal[separate],xerr=False,emptybins=False,axes=ax,label="Jet axis reference",fmt='o') #Plot jT histogram, 
    rplt.errorbar(signal2[separate],xerr=False,emptybins=False,axes=ax,label="Leading track reference",fmt='o')
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(jetPt[separate][0],jetPt[separate][1])) 
    ax.set_xlim([0.1,12])
    ax.set_ylim([5e-6,1.5e3])
    ax.legend(loc = 'lower left')
    
    plt.savefig("PythonFigures/MixedFullJetsR04JetConeJtLeadingRefJetPt{0}.pdf".format(separate),format='pdf') #Save figure
    plt.show() #Draw figure on screen

  else:
    n_rows = n_figs//4
    print(n_rows)
    fig, axs = defs.makegrid(4,2,xlog=True,ylog=True,d=d,shareY=False,figsize=(10,5))
    axs = axs.reshape(8)
    #axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone',fontsize = 7)
    ratios = []
    for jT,jT2,pT,ax,i in zip(JtJet[start:],JtLeading[start:],jetPt[start:],axs[0:4],range(0,9)):
      jT.SetMarkerColor(1)
      jT.SetMarkerStyle(24)
      jT.SetLineColor(1)
      #jT.SetMarkerSize(mSize)
      jT2.SetMarkerColor(2)
      jT2.SetMarkerStyle(25)
      #jT2.SetMarkerSize(mSize)
      jT2.SetLineColor(2)
      ratio = jT2.Clone()
      ratio.Divide(jT)
      ratios.append(ratio)
      plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label="Jet axis reference",fmt='o',fillstyle='none',ecolor='black') #Plot jT histogram, 
      line = plot.get_children()[0]
      line.set_markersize(mSize)
      line.set_markerfacecolor('none')

      plot = rplt.errorbar(jT2,xerr=False,emptybins=False,axes=ax,label="Leading track reference",fmt='o',fillstyle='none',ecolor='red') #Plot jT histogram, 
      line = plot.get_children()[0]
      line.set_markersize(mSize)
      line.set_markerfacecolor('none')
      #line.set_color(color)
      
      ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
  
      ax.set_xlim([0.1,22]) #Set x-axis limits
      ax.set_ylim([5e-5,2e3]) #Set y-axis limits
      ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
    
    for ratio,ax,i in zip(ratios,axs[4:],range(0,9)):
      plot = rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,fmt='o',fillstyle='none',ecolor='black') #Plot jT histogram, 
      line = plot.get_children()[0]
      line.set_markersize(mSize)
      line.set_markerfacecolor('none')
      ax.set_yscale('linear')
      ax.set_xlim([0.1,22]) #Set x-axis limits
      ax.set_ylim([0,5]) #Set y-axis limits
    
    handles, labels = axs[0].get_legend_handles_labels()
    handles = [container.ErrorbarContainer(h,has_xerr=False,has_yerr=True) if isinstance(h, container.ErrorbarContainer) else h for h in handles]
    axs[0].legend(handles,labels,loc = 'lower left',numpoints=1)
    axs[4].set_ylabel('Ratio')
    axs[7].set_ylabel('Ratio')
        
    plt.savefig("PythonFigures/MixedFullJetsR04JetConeJtLeadingRefPtFrom{}To{}.pdf".format(start,end),format='pdf') #Save figure
    plt.show() #Draw figure on screen


def drawWithErrors2Combined(hists,hists2,titles,xlow,xhigh,logx,ylow,yhigh,ylog,xTitle,yTitle,title,file):
  colors = ['white','black','red','green','blue']

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
  ax = axs[0]
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  ax.set_xlabel(xTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  ax.set_ylabel(yTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  axs[1].set_ylabel(yTitle,fontsize=labelsize)
  ax.set_ylabel(yTitle,fontsize=labelsize) #Add x-axis labels for bottom row

  ratios1 = []
  ratios2 = []
  if logx:
    ax.set_xscale('log')
  for h,h2,c,R in zip(hists,hists2,(1,2,3,4),titles):
    ax = axs[0]
    h.SetMarkerColor(c)
    h.SetMarkerStyle(24)
    h2.SetMarkerColor(c)
    h2.SetMarkerStyle(25)
    ratios1.append(grrDivide(h, hists[0]))
    ratios2.append(grrDivide(h2, hists2[0]))
    plot = rplt.errorbar(h,xerr=False,emptybins=False,label = R,axes=ax,fmt='+')

    line = plot.get_children()[0]
    ax.set_xlim([xlow,xhigh])
    ax.set_ylim([ylow,yhigh])
    ax = axs[1]
    ax.grid(True)
    ax.set_xlim([xlow,xhigh])
    ax.set_ylim([ylow,yhigh])
    plot = rplt.errorbar(h2,xerr=False,yerr=True,emptybins=False,axes=ax,label=R,fmt='+')

  
  print("({} - {})/9.0 + {} is {}".format(yhigh,ylow,ylow,(yhigh-ylow)/9.0+ylow))
  axs[1].text(65,0.2,title + '\n' + d['system'] +'\n'+  d['jettype'] + '\n' + r'Anti-$k_T$',fontsize = 10)
  axs[1].text(102,1.22,"Wide")
  axs[0].text(102,1.22,"Narrow")


  axs[0].legend(loc = 'upper left')
  for ax in axs[0:2]:
    ax.set_xlim([xlow,xhigh])
    ax.set_ylim([ylow,yhigh])
    ax.grid(True)

  
  ax = axs[2]
  ax2 = axs[3]
  ax.grid(True)
  ax2.grid(True)
  ax.set_xlabel(xTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  ax2.set_xlabel(xTitle,fontsize=labelsize) #Add x-axis labels for bottom row
  ax.set_ylabel('Ratio',fontsize=labelsize) #Add x-axis labels for bottom row
  ax2.set_ylabel('Ratio',fontsize=labelsize) #Add x-axis labels for bottom row

  for ratio,ratio2,c in zip(ratios1,ratios2,(1,2,3,4)):
    ratio.SetMarkerColor(c)
    ratio2.SetMarkerColor(c)
    ratio.SetMarkerStyle(24)
    ratio2.SetMarkerStyle(25)
    plot = rplt.errorbar(ratio,xerr = False,yerr=True,emptybins=False,axes=ax)
    plot = rplt.errorbar(ratio2,xerr = False,yerr=True,emptybins=False,axes=ax2)
        
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([0.1,2.6])
  ax2.set_xlim([xlow,xhigh])
  ax2.set_ylim([0.1,2.6]) 
  ax2.grid(True) 
  if logx:
    ax.set_xscale('log')
    ax2.set_xscale('log')
  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0  
  plt.savefig("{}.pdf".format(file),format='pdf') #Save figure

  plt.show() #Draw figure on screen
   
 
  
if __name__ == "__main__": main()
