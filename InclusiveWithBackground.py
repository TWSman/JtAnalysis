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
  print 'Number of arguments: ', len(sys.argv), 'arguments.'
  print 'Argument list:',str(sys.argv)
  filename = sys.argv[1]
  print "Input file: "
  print filename
  #FullJets_R04 = dataset("FullR04",NFIN=0,filename=filename,directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=5)
  Mixed_FullJets_R04 = datasetMixed("FullR04",NFIN=0,range=5,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=2)
  JtWithBackgroundRatio(Mixed_FullJets_R04, 'JetConeJtWeightBin', 'BgJtWeightBin','BgRndmJtWeightBin',includeRandom=True,start=4,step=1)
  plt.savefig("PythonFigures/MixedFullJetsR04JetConeJtInclusiveComparison.pdf",format='pdf') #Save figure
  plt.show() #Draw figure on screen


  
if __name__ == "__main__": main()
