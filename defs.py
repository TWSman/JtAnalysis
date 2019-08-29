import rootpy
import re
import matplotlib
from rootpy.plotting import Canvas,Graph,Legend,Hist
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import LogLocator
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
from ROOT import TF1,TMath,Double
import numpy as np
import math
labelsize= 15
from matplotlib.gridspec import GridSpec


B1low = 0.1
B1start = [0.15,0.15,0.17,0.18,0.15,0.15,0.21,0.15,0.15]
B1high = 0.5
B2low = 40.0
B2high = 150.0
B2start = [50,50,66.20,67.90,50,50,75.89,77.62,77.62]
B3cut = 6
B3low = 0.0
B3high1 = 10.0
B3high2 = 25.0
B3start = [7,7,5.06,4.90,10,10,7,12.93,12.93]
B4low = 2.0
B4high = 15.0
B4start = [3,3,9.91,8.80,3,3,5.62,4.18,4.18]
B5low = 0.90
B5high = 4.5
B5start = [1.5,1.5,4.5,4.5,1.5,1.5,2.88,1.62,1.62]
peakstart = 0.5
peaklow = 0.3
peakhigh = 0.8
Rs = (0.3,0.4,0.5,-0.4)

def makegrid(nx=4,ny=2,xlog=True,ylog=True,d=None,shareY=True,figsize = (10,5),**kwargs):
  """Create and nx by ny grid of subfigures with shared y-axes
  Args:
    nx: Number of figures in y direction, default is 4
    ny: Number of figures in x direction, deftault is 2
    xlog: Whether to use logarithmic scale on X-axis, default is True
    ylog: Whether to use logarithmic scale on Y-axis, default is True
    d: Not used
    shareY: Whether subplots should share Y-axis, default is True
  """
  if('xtitle' in kwargs):
    xtitle = kwargs['xtitle']
  else:
    xtitle = r'$j_\mathrm{T}\left(\mathrm{GeV}/c\right)$'
  if('ytitle' in kwargs):
    ytitle = kwargs['ytitle']
  else:
    ytitle = r'$\frac{1}{N_{jets}}\frac{\mathrm{d} N}{j_\mathrm{T}\mathrm{d} j_\mathrm{T}}$'

  fig, axs = plt.subplots(ny,nx,figsize=figsize,sharey=shareY,sharex=True) #Create figure with 8 subfigures, axs is a list of subfigures, fig is the whole thing
  #axs = axs.reshape(nx*ny) #Because the figures is in a nx*ny layout axs is a 2 dimensional array with nx * ny elements, this makes it a 1 dimensional array with nx*ny  elements
  #axs[0][0].text(0.02,0.005,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n Charged jT\n' r'Anti-$k_T$, R=0.4' '\nJet Cone',fontsize=7) #Add text to second subfigure, first parameters are coordinates in the drawn scale/units
  #axs[1][0].text(0.02,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'],fontsize = 7)
#   for ax, i in zip(axs,range(nx*ny)):
#     if(i%nx == 0 or i % nx == nx-1):
#       ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{\mathrm{d} N}{j_\mathrm{T}\mathrm{d} j_\mathrm{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
#     if(i % nx == nx-1):
#       ax.yaxis.set_label_position('right') #Set the y-axis label position to right hand side for the rightmost subfigures
#   for ax in axs[-nx:]:
#     ax.set_xlabel(r'$j_\mathrm{T}\left[GeV\right]$') #Add x-axis labels for bottom row
  if(ny == 1):
    axs0 = [axs]
  else:
    axs0 = axs
  for ax in axs0:
    ax[0].set_ylabel(ytitle,fontsize=18) #Add y-axis labels to left- and righmost subfigures
    if(nx > 2):
      ax[-1].set_ylabel(ytitle,fontsize=18) #Add y-axis labels to left- and righmost subfigures
      ax[-1].yaxis.set_label_position('right') #Set the y-axis label position to right hand side for the rightmost subfigures
    for a in ax[1:]:
      a.yaxis.tick_right()  
  
  for ax in axs0[-1]:
    ax.set_xlabel(xtitle,fontsize=15) #Add x-axis labels for bottom row

  for axs1 in axs0: 
    for ax in axs1:
      if(xlog):
        ax.set_xscale('log') #Set logarithmic scale
      if(ylog):
        ax.set_yscale('log')
      ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
      ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
      ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
      #ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
      ax.set_xlim([0.1,20]) #Set x-axis limits
      ax.set_ylim([5e-4,2e3]) #Set y-axis limits
      #plt.setp(ax.get_xticklabels()[-3],visible = False) #Remove last label
      ax.grid(True) #Draw grid
      
  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
  #print axs.ndim
  return fig,axs


  
def makeRatio(xlog=True,ylog=True,d=None,shareY=False,figsize = (4,6),**kwargs):
  """
  
  Args:

    xlog: Whether to use logarithmic scale on X-axis, default is True
    ylog: Whether to use logarithmic scale on Y-axis, default is True
    d: Not used
    shareY: Whether subplots should share Y-axis, default is False
  """
  if('xtitle' in kwargs):
    xtitle = kwargs['xtitle']
  else:
    xtitle = r'$j_\mathrm{T}\left(\mathrm{GeV}/c\right)$'
  if('ytitle' in kwargs):
    ytitle = kwargs['ytitle']
  else:
    ytitle = r'$\frac{1}{N_{jets}}\;\frac{1}{j_\mathrm{T}}\;\frac{\mathrm{d} N}{\mathrm{d} j_\mathrm{T}}$'

  fig=plt.figure(figsize=figsize,)
  gs=GridSpec(2,1,height_ratios=[2, 1])
  ax0 = plt.subplot(gs[0])
  ax1 = plt.subplot(gs[1])
  axs = np.array([ax0,ax1])
  
  #fig, axs = plt.subplots(2,1,figsize=figsize,sharey=shareY,sharex=True) 
  #axs = axs.reshape(nx*ny) #Because the figures is in a nx*ny layout axs is a 2 dimensional array with nx * ny elements, this makes it a 1 dimensional array with nx*ny  elements
  #axs[0][0].text(0.02,0.005,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n Charged jT\n' r'Anti-$k_T$, R=0.4' '\nJet Cone',fontsize=7) #Add text to second subfigure, first parameters are coordinates in the drawn scale/units
  #axs[1][0].text(0.02,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'],fontsize = 7)
#   for ax, i in zip(axs,range(nx*ny)):
#     if(i%nx == 0 or i % nx == nx-1):
#       ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{\mathrm{d} N}{j_\mathrm{T}\mathrm{d} j_\mathrm{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
#     if(i % nx == nx-1):
#       ax.yaxis.set_label_position('right') #Set the y-axis label position to right hand side for the rightmost subfigures
#   for ax in axs[-nx:]:
#     ax.set_xlabel(r'$j_\mathrm{T}\left[GeV\right]$') #Add x-axis labels for bottom row
  axs[0].set_ylabel(ytitle,fontsize=16) #Add y-axis labels to left- and righmost subfigures

  for ax in axs:
    ax.set_xlabel(xtitle,fontsize=15) #Add x-axis labels for bottom row
    if(xlog):
      ax.set_xscale('log') #Set logarithmic scale
    if(ylog):
      ax.set_yscale('log')

    ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
    ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
    ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
    #ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
    ax.set_xlim([0.1,20]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    #plt.setp(ax.get_xticklabels()[-3],visible = False) #Remove last label
    if(kwargs.get('grid',True)):
      ax.grid(True) #Draw grid      
  axs[0].xaxis.set_ticklabels([])
  #axs[0].xaxis.set_ticks([])
  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
  #print axs.ndim
  return fig,axs



def make8grid(xlog=True,ylog=True,bins=None,d=None):
  """Create an 4 by 2 grid of subfigures with shared axes
  
  Args:
    xlog: Whether to use logarithmic scale on X-axis, default is True
    ylog: Whether to use logarithmic scale on Y-axis, default is True
    d: Not used
    bins: Not used
  """
  fig, axs = plt.subplots(2,4,figsize=(10,5),sharey=True,sharex=True) #Create figure with 8 subfigures, axs is a list of subfigures, fig is the whole thing
  axs = axs.reshape(8) #Because the figures is in a 2x4 layout axs is a 2 dimensional array with 2x4 elements, this makes it a 1 dimensional array with 8 elements

  axs[1].text(0.02,0.005,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n Charged jT\n' r'Anti-$k_T$, R=0.4' '\nJet Cone',fontsize=7) #Add text to second subfigure, first parameters are coordinates in the drawn scale/units
  axs[2].text(0.02,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'],fontsize = 7)
  for ax in [axs[0],axs[3],axs[4],axs[7]]: 
    ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{\mathrm{d} N}{j_\mathrm{T}\mathrm{d} j_\mathrm{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
  for ax in axs[4:]:
    ax.set_xlabel(r'$j_\mathrm{T}\left(\mathrm{GeV}/c\right)$',fontsize=15) #Add x-axis labels for bottom row
  for ax in [axs[3],axs[7]]: 
    ax.yaxis.set_label_position('right') #Set the y-axis label position to right hand side for the rightmost subfigures

  for (ax,pT) in zip (axs,bins): 
    if(xlog):
      ax.set_xscale('log') #Set logarithmic scale
    if(ylog):
      ax.set_yscale('log')
    ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
    ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
    ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
    ax.set_xlim([0.1,20]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    plt.setp(ax.get_xticklabels()[-3],visible = False) #Remove last label
    ax.grid(True) #Draw grid 

  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
  return fig,axs

#def fitJtHisto(measJt):
  


def draw8grid(measJt,measBg,jetPt,xlog = True,ylog = True,name="newfile.pdf"):
  """Create an 4 by 2 grid of subfigures with shared axes and plots jT with background in jet pT bins
  
  Args:
    measJt: List of jT histograms
    measBg: List of Bg Histograms
    jetPt: List of jet Pt bins in tuples (low border, high border)
    xlog: Whether to use logarithmic scale on X-axis, default is True
    ylog: Whether to use logarithmic scale on Y-axis, default is True
    name: Name of output file
  """
  fig, axs = plt.subplots(2,4,figsize=(10,5),sharey=True,sharex=True) #Create figure with 8 subfigures, axs is a list of subfigures, fig is the whole thing
  axs = axs.reshape(8) #Because the figures is in a 2x4 layout axs is a 2 dimensional array with 2x4 elements, this makes it a 1 dimensional array with 8 elements

  axs[1].text(0.02,0.005,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n Charged jT\n' r'Anti-$k_T$, R=0.4' '\nJet Cone',fontsize=7) #Add text to second subfigure, first parameters are coordinates in the drawn scale/units
  for ax in [axs[0],axs[3],axs[4],axs[7]]: 
    ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{\mathrm{d} N}{j_\mathrm{T}\mathrm{d} j_\mathrm{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
  for ax in axs[4:]:
    ax.set_xlabel(r'$j_\mathrm{T}\left(\mathrm{GeV}/c\right)$',fontsize=15) #Add x-axis labels for bottom row
  for ax in [axs[3],axs[7]]: 
    ax.yaxis.set_label_position('right') #Set the y-axis label position to right hand side for the rightmost subfigures

  for (jT,Bg,ax,i,pT) in zip (measJt,measBg,axs,range(0,8),jetPt): 
    if(xlog):
      ax.set_xscale('log') #Set logarithmic scale
    if(ylog):
      ax.set_yscale('log')
    rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label='jT',fmt='+') #Plot jT histogram, 
    rplt.errorbar(Bg,xerr=False,emptybins=False,axes=ax,label='jT Bg',fmt='go') #Plot bg jT histogram
    if i == 0: #For the first subfigure add a legend to bottom left corner
      ax.legend(loc ='lower left')
    ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
    ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
    ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
    ax.set_xlim([0.1,5]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    ax.grid(True) #Draw grid 

  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
  plt.savefig(name,format='pdf') #Save figure
  plt.show() #Draw figure on screen
  
def fitJtHisto(histo,method,cut,iJet,iFinder,title="",draw=False):
  d = {}
  gaussfit = TF1("gaussfit","gausn",0,10)
  gaussfit.FixParameter(1,0)
  print("Set parameter B1 limits to {}-{}".format(B1low,B1high))
  gaussfit.SetParLimits(2,B1low,B1high)
  print("Set parameter B2 to {}".format(B2start[iJet]))
  gaussfit.SetParameter(0,B2start[iJet])
  print("Set parameter B2 limits to {}-{}".format(B2low,B2high))
  gaussfit.SetParLimits(0,B2low,B2high)
  histo.Fit("gaussfit","QN")
  if method == "alt":
    invG = TF1("invG","[0]* (pow([1]*([2]+1),[2])/TMath::Gamma([2]))*exp(-[1]*([2]+1)/x) * pow(x, -[2]-1)",0,10)
  else:
    invG = TF1("invG","[0]* (pow([1],[2])/TMath::Gamma([2]))*exp(-[1]/x) * pow(x, -[2]-1)",0,10)
  print("Set parameter B3 start to {}".format(B3start[iJet]))
  invG.SetParameter(0,B3start[iJet])
  if iJet < B3cut:
    print("Set parameter B3 limits to {}-{}".format(B3low,B3high1))
    invG.SetParLimits(0,B3low,B3high1)
  else:
    print("Set parameter B3 limits to {}-{}".format(B3low,B3high2))
    invG.SetParLimits(0,B3low,B3high2)
  if method == "alt":
    invG.SetParameter(1,peakstart)
    invG.SetParLimits(1,peaklow,peakhigh)
  else:
    print("Set parameter B5 start to {}".format(B5start[iJet]))
    invG.SetParameter(1,B5start[iJet])
    print("Set parameter B5 limits to {}-{}".format(B5low,B5high))
    invG.SetParLimits(1,B5low,B5high)
  
  print("Set parameter B4 start to {}".format(B4start[iJet]))
  invG.SetParameter(2,B4start[iJet])
  print("Set parameter B4 limits to {}-{}".format(B4low,B4high))
  invG.SetParLimits(2,B4low,B4high)
  histo.Fit("invG","QN","",cut,3)
  if method == "alt":
    gaussfit3= TF1("gaussfit3","gausn(0) + [3]* (pow([4]*([5]+1),[5])/TMath::Gamma([5]))*exp(-[4]*([5]+1)/x) * pow(x, -[5]-1)",0,10)
  else:
    gaussfit3= TF1("gaussfit3","gausn(0) + [3]* (pow([4],[5])/TMath::Gamma([5]))*exp(-[4]/x) * pow(x, -[5]-1)",0,10)
  gaussfit3.SetParameter(0, gaussfit.GetParameter(0))
  gaussfit3.SetParameter(1, gaussfit.GetParameter(1))
  gaussfit3.SetParameter(2, gaussfit.GetParameter(2))
  gaussfit3.SetParameter(3, invG.GetParameter(0))
  gaussfit3.SetParameter(4, invG.GetParameter(1))
  gaussfit3.SetParameter(5, invG.GetParameter(2))
  gaussfit3.FixParameter(1,0)

  gaussfit3.SetParLimits(0,B2low,B2high) #B2
  gaussfit3.SetParLimits(2,B1low,B1high) #B1
  if iJet < B3cut:
    gaussfit3.SetParLimits(3,B3low,B3high1) #B3
  else:
    gaussfit3.SetParLimits(3,B3low,B3high2) #B3
  if method == "alt":
    gaussfit3.SetParLimits(4,peaklow,peakhigh) #Peak
  else:
    gaussfit3.SetParLimits(4,B5low,B5high) #B5
  gaussfit3.SetParLimits(5,B4low,B4high) #B4

  lastBin = histo.FindLastBinAbove(1e-5)
  end = histo.GetBinCenter(lastBin)
  histo.Fit("gaussfit3","QN","",0.1,end)
  #chi2[iF][ij] = gaussfit3.GetChisquare()
  #chi2[iF][ij] = getChi2(histo,gaussfit3,0.1,end,0)
  #chi2dof[iF][ij] = chi2[iF][ij]/(bins[iF][ij]-5)
  #chi2n[iF][ij] = chi2[iF][ij]/bins[iF][ij]
  B2 = gaussfit3.GetParameter(0)
  B1 = gaussfit3.GetParameter(2)
  B2e = gaussfit3.GetParError(0)
  B1e = gaussfit3.GetParError(2)
  gaussigma = TMath.Sqrt(2)*B1
  gausyield = B2*B1/TMath.Sqrt(2 * TMath.Pi())
  gausyielde = TMath.Sqrt((B2*B2*B1e*B1e + B1*B1*B2e*B2e)/(2*TMath.Pi()))
  gaussigmae = TMath.Sqrt(2)*B1e

  constant = gaussfit3.GetParameter(3) #B3
  constante = gaussfit3.GetParError(3)
  alpha = gaussfit3.GetParameter(5) # B4
  alphae = gaussfit3.GetParError(5)
  if method == "alt":
    peak = gaussfit3.GetParameter(4)
    peake = gaussfit3.GetParError(4)
    beta = peak*(alpha+1)#B5
    betae = TMath.Sqrt(pow( (alpha+1) * peake,2) + pow(peak*alphae,2))
  else:
    beta = gaussfit3.GetParameter(4) # B5
    betae = gaussfit3.GetParError(4)
    peak = beta/(alpha+1)
    peake = TMath.Sqrt(pow(betae/(alpha+1),2)+pow(alphae*beta/pow(alpha+1,2),2))
    
  gammaYield = constant * beta /(alpha - 1)
  gammaRMS = beta /TMath.Sqrt((alpha-2)*(alpha-3))
  gammaYielde =TMath.Sqrt(pow(beta *  constante/(alpha-1),2) + pow(constant*beta*alphae/pow(alpha-1,2),2) + pow(constant * betae/(alpha-1),2))
  gammaRMSe = TMath.Sqrt(pow((5-2*alpha)*beta*alphae/pow(2*((alpha-2)*(alpha-3)),1.5),2) + pow(betae/TMath.Sqrt((alpha-2)*(alpha-3)),2))

  d['B2'] = B2
  d['B1'] = B1
  d['B3'] = constant
  d['B4'] = alpha
  d['B5'] = beta
  d['peak'] = peak
  d['B2e'] = B2e
  d['B1e'] = B1e
  d['B3e'] = constante
  d['B4e'] = alphae
  d['B5e'] = betae
  d['peake'] = peake

  d['gausRMS'] = gaussigma
  d['gausRMSe'] = gaussigmae
  d['gausYield'] = gausyield
  d['gausYielde'] = gausyielde
  d['gammaRMS'] = gammaRMS
  d['gammaRMSe'] = gammaRMSe
  d['gammaYield'] = gammaYield
  d['gammaYielde'] = gammaYielde
  if(draw):
    drawFit(histo,gaussfit3,iJet,iFinder,title)
  return gaussfit3,d

def drawFit(histo,fit,iJet,iFinder,title):
  jetPt = [5,10,20,30,40,60,80,100,150,500]
  fig = plt.figure(1)
  ax = fig.add_subplot(1,1,1)
  ax.set_xscale('log')
  ax.set_yscale('log')
  #ax.text(0.2,0.0005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone',fontsize = 7)
  rplt.errorbar(histo,xerr=False,emptybins=False,axes=ax,label="jT",fmt='o') #Plot jT histogram,
  xs = np.arange(0,10,0.01).tolist()
  #print(xs[0:50])
  #print(np.shape(xs))
  ys = [fit.Eval(x) for x in xs]
  #print(type(xs))
  #print(type(ys))
  #print(ys[0:50])
  #print(np.shape(ys))
  #print(ys[0][0:50])
  #print(type(ys[0]))
  #print("moi")
  #rplt.hist(Hist(fit.GetHistogram()),axes=ax,label = "Fit",fmt='-')
  plt.plot(xs,ys)
  ax.text(0.3,1e2,title + '\n' + r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(jetPt[iJet],jetPt[iJet+1])) 
  ax.set_xlim([0.1,12])
  ax.set_ylim([5e-6,1.5e3])
  ax.legend(loc = 'lower left')
  
  #plt.savefig("PythonFigures/UnfoldedFullJetsR04JetConeJtFinalJetPt{0}.pdf".format(separate),format='pdf') #Save figure
  plt.show() #Draw figure on screen

def grrScale(gr1,scale):
  x1 = Double()
  x_ = Double()
  x1e = Double()

  y1 = Double()
  y1e = Double()

  y = []
  ex = []
  ey = []

  NC =  gr1.GetN()
  x = [0 for i in range(NC)]

  for ii in range(NC):
    x1e = gr1.GetErrorX(ii)
    y1e = gr1.GetErrorY(ii)
    gr1.GetPoint(ii, x_, y1)
    x1 = x_ * 1.0
    x[ii] = x1
    y.append(y1*scale)
    ex.append(x1e)
    ey.append(y1e*scale)
    
  gr= Graph(NC)
  for x0,y0,x0e,y0e,i in zip(x,y,ex,ey,range(NC)):
    gr.SetPoint(i,x0,y0)
    gr.SetPointError(i,x0e,x0e,y0e,y0e)
  return gr


def grrDivide(gr1,gr2):
  x1 = Double()
  x_ = Double()
  x1e = Double()
  x2 = Double()
  x2e = Double()
  y1 = Double()
  y1e = Double()
  y2 = Double()
  y2e = Double()
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


def makeSystError(gr1,gr2,**kwargs):
  x_ = Double()
  x1e = Double()
  x2 = Double()
  x2e = Double()
  y1 = Double()
  y1e = Double()
  y2 = Double()
  y2e = Double()
  NC = gr1.GetN()
  print("Number of points: {}".format(NC))
  y = []
  ex = []
  ey = []
  x = [0 for i in range(NC)]
  
  for ii in range(NC):
    x1e = gr1.GetErrorX(ii)
    y1e = gr1.GetErrorY(ii)
    gr1.GetPoint(ii, x_, y1)
    x1 = x_ * 1.0
    gr2.GetPoint(ii, x2, y2)
    x2e = gr2.GetErrorX(ii)
    y2e = gr2.GetErrorY(ii)
    x[ii] = x1
    if(kwargs.get('abs',False)):
      y.append(math.fabs(y1-y2))
    elif(kwargs.get('rel',False)):
      y.append((y1-y2)/(y1+y2))
    else:
      y.append(y1-y2)
    ex.append(x1e)
    ey.append(math.sqrt(y1e**2+y2e**2))
  gr= Graph(NC)
  for x0,y0,x0e,y0e,i in zip(x,y,ex,ey,range(NC)):
    gr.SetPoint(i,x0,y0)
    print("Set point {} {},{}".format(i,x0,y0))
    gr.SetPointError(i,x0e,x0e,y0e,y0e)
  return gr

def addConstantError(gr1,error,**kwargs):
  x_ = Double()
  x1e = Double()
  y1_ = Double()
  y1e = Double()
  NC = gr1.GetN()
  y2 = [0 for i in range(NC)]
  ex = []
  ey = []
  x = [0 for i in range(NC)]
  gr= Graph(NC)
  for ii in range(NC):
    x1e = gr1.GetErrorX(ii)
    y1e = gr1.GetErrorY(ii)
    gr1.GetPoint(ii, x_, y1_)
    print("{} {} {}".format(i,x_,y1_))
    x1 = x_ * 1.0
    
    x[ii] = x1
    y1 = y1_ * 1.0
    y2[ii] = y1
    print("Append y with {}".format(y1))
    print(y2)
    ex.append(x1e)
    if(kwargs.get('rel'),False):
      ey.append(math.sqrt(y1e**2+(error*y1)**2))
      print("Old error: {}, New Error: {}".format(y1e,math.sqrt(y1e**2+(error*y1)**2)))
    else:
      ey.append(math.sqrt(y1e**2+error**2))
#     gr.SetPoint(i,x1,y1)
#     print("Set point {} {},{}".format(i,x1,y1))
# 
#     if(kwargs.get('rel'),False):
#       gr.SetPointError(i,x1e,x1e,math.sqrt(y1e**2+(error*y1)**2),math.sqrt(y1e**2+(error*y1)**2))
#       print("Set point Error {} {},{}".format(i,x1e,math.sqrt(y1e**2+(error*y1)**2)))
#     else:
#       gr.SetPointError(i,x1e,x1e,math.sqrt(y1e**2+error**2),math.sqrt(y1e**2+error**2))

  gr= Graph(NC)
  print(x)
  print(y2)
  print(ex)
  print(ey)
  for x0,y0,x0e,y0e,i in zip(x,y2,ex,ey,range(NC)):
    gr.SetPoint(i,x0,y0)
    print("Set point {} {},{}".format(i,x0,y0))
    gr.SetPointError(i,x0e,x0e,y0e,y0e)
    print("Set point Error {} {},{}".format(i,x0e,y0e))

  print(gr)
  return gr


def drawErrors2(error,xlow,  xhigh, xlog, ylow, yhigh, ylog, title,title2,comment,file,iS, ij,limit,error2=None,title3=""):
  jetPt = [5,10,20,30,40,60,80,100,150,500]
  fig = plt.figure(1)
  ax = fig.add_subplot(1,1,1)
  if(xlog):
    ax.set_xscale('log')
  if(ylog):
    ax.set_yscale('log')
  ax.set_xlabel(r'$p_{{T,\mathrm{{jet}}}}$',fontsize=labelsize) #Add x-axis labels for bottom row
  ax.set_ylabel("Relative error",fontsize=labelsize)
  #ax.text(0.2,0.0005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone',fontsize = 7)
  rplt.errorbar(error,xerr=False,emptybins=False,axes=ax,label=title2,fmt='o') #Plot jT histogram,
  if(error2 is not None):
    error2.SetMarkerColor(2)
    rplt.errorbar(error2,xerr=False,emptybins=False,axes=ax,label=title3,fmt='o') #Plot jT histogram,    
  xs = np.arange(xlow,xhigh,0.01).tolist()
  ys = [limit for x in xs]
  ys2 = [-1 * limit for x in xs]

  #rplt.hist(Hist(fit.GetHistogram()),axes=ax,label = "Fit",fmt='-')
  plt.plot(xs,ys,linestyle='dashed',color='blue',label='{}%'.format(limit*100))
  plt.plot(xs,ys2,linestyle='dashed',color='blue')

  if(ij > 0):
    ax.text(xhigh*0.8,3*(yhigh+ylow)/4,title + '\n' + r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(jetPt[iJet],jetPt[iJet+1]))
  else:
    ax.text(xlow+(xhigh-xlow)*0.1,0.8*yhigh,title)
  ax.set_xlim([xlow,xhigh])
  ax.set_ylim([ylow,yhigh])
  ax.legend(loc = 'upper right')
  print("Save file {}".format(file))
  plt.savefig(file,format='pdf') #Save figure
  plt.show() #Draw figure on screen
