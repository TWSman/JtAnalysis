import rootpy
import re
import matplotlib
from rootpy.plotting import Canvas
from rootpy.plotting import Legend
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import LogLocator
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt

def makegrid(nx=4,ny=2,xlog=True,ylog=True,d=None,shareY=True,figsize = (10,5)):
  """Create and nx by ny grid of subfigures with shared y-axes
  
  Args:
    nx: Number of figures in y direction, default is 4
    ny: Number of figures in x direction, deftault is 2
    xlog: Whether to use logarithmic scale on X-axis, default is True
    ylog: Whether to use logarithmic scale on Y-axis, default is True
    d: Not used
    shareY: Whether subplots should share Y-axis, default is True
  """
  fig, axs = plt.subplots(ny,nx,figsize=figsize,sharey=shareY,sharex=True) #Create figure with 8 subfigures, axs is a list of subfigures, fig is the whole thing
  #axs = axs.reshape(nx*ny) #Because the figures is in a nx*ny layout axs is a 2 dimensional array with nx * ny elements, this makes it a 1 dimensional array with nx*ny  elements
  #axs[0][0].text(0.02,0.005,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n Charged jT\n' r'Anti-$k_T$, R=0.4' '\nJet Cone',fontsize=7) #Add text to second subfigure, first parameters are coordinates in the drawn scale/units
  #axs[1][0].text(0.02,0.005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'],fontsize = 7)
#   for ax, i in zip(axs,range(nx*ny)):
#     if(i%nx == 0 or i % nx == nx-1):
#       ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{dN}{j_{T}dj_{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
#     if(i % nx == nx-1):
#       ax.yaxis.set_label_position('right') #Set the y-axis label position to right hand side for the rightmost subfigures
#   for ax in axs[-nx:]:
#     ax.set_xlabel(r'$j_{T}\left[GeV\right]$') #Add x-axis labels for bottom row
  for ax in axs:
    ax[0].set_ylabel(r'$\frac{1}{N_{jets}}\frac{dN}{j_{T}dj_{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
    ax[-1].set_ylabel(r'$\frac{1}{N_{jets}}\frac{dN}{j_{T}dj_{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
    ax[-1].yaxis.set_label_position('right') #Set the y-axis label position to right hand side for the rightmost subfigures
    for a in ax[1:]:
      a.yaxis.tick_right()  
  
  for ax in axs[-1]:
    ax.set_xlabel(r'$j_{T}\left[GeV\right]$') #Add x-axis labels for bottom row

  for axs1 in axs: 
    for ax in axs1:
      if(xlog):
        ax.set_xscale('log') #Set logarithmic scale
      if(ylog):
        ax.set_yscale('log')
      ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
      ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
      ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
      #ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
      ax.set_xlim([0.01,20]) #Set x-axis limits
      ax.set_ylim([5e-4,2e3]) #Set y-axis limits
      plt.setp(ax.get_xticklabels()[-3],visible = False) #Remove last label
      ax.grid(True) #Draw grid
      
  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
  print axs.ndim
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
    ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{dN}{j_{T}dj_{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
  for ax in axs[4:]:
    ax.set_xlabel(r'$j_{T}\left[GeV\right]$') #Add x-axis labels for bottom row
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
    ax.set_xlim([0.01,20]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    plt.setp(ax.get_xticklabels()[-3],visible = False) #Remove last label
    ax.grid(True) #Draw grid 

  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
  return fig,axs

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
    ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{dN}{j_{T}dj_{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
  for ax in axs[4:]:
    ax.set_xlabel(r'$j_{T}\left[GeV\right]$') #Add x-axis labels for bottom row
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
