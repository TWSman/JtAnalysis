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
  MB_FullJets_R04 = dataset("MBFullR04",NFIN=0,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',color=1,style=24,rebin=5)
  Triggered_FullJets_R04 = dataset("TriggeredFullR04",NFIN=0,filename=filename,directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=5)
  compareSetsWithRatio((MB_FullJets_R04,Triggered_FullJets_R04),'JetConeJtWeightBin')  
  plt.savefig("PythonFigures/MBvsTriggeredFullJetsR04JetConeJt.pdf",format='pdf') #Save figure
  plt.show() #Draw figure on screen


  
if __name__ == "__main__": main()
