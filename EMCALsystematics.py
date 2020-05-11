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
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
from dataset import *
def main(): 
  print('Number of arguments: ', len(sys.argv), 'arguments.')
  print('Argument list:',str(sys.argv))
  filename1 = "Pythia/pythia8226_pp5TeV_Monash2013_clusterCorr100_9000.root"
  filename2 = "Pythia/pythia8226_pp5TeV_Monash2013_clusterCorr98_6800.root"
  filename3 = "Pythia/pythia8226_pp5TeV_Monash2013_clusterCorr102_5000.root"
  filename4 = "Pythia/pythia8226_pp5TeV_Monash2013_JetScaleUp_4000.root"
  filename5 = "Pythia/pythia8226_pp5TeV_Monash2013_JetScaleDown_2000.root"
  filename_default = "JetJtAna_Default_20190121_153453.root"
  filename_hadcorr = "JetJtAna_HadCorr_20190121_152726.root"

#   if(len(sys.argv) > 3):
#     start = int(sys.argv[3])
#   else:
#     start = 1
  start = 2

  Default_data = dataset("LHC13b_default",NFIN=0,range=(0,6),filename=filename_default,directory='AliJJetJtTask/AliJJetJtHistManager',color=2,style=24,rebin=4)
  HadCorr_data = dataset("LHC13b_hadcorr",NFIN=0,range=(0,6),filename=filename_hadcorr,directory='AliJJetJtTask/AliJJetJtHistManager',color=2,style=24,rebin=4)

  Default = dataset("Default",NFIN=0,range=(0,8),filename=filename1,directory='JCDijetBaseTask/jcdijet',color=2,style=24,rebin=4)
  HadCorr_high = dataset("-2%",NFIN=0,range=(0,8),filename=filename2,directory='JCDijetBaseTask/jcdijet',color=2,style=24,rebin=4)
  HadCorr_low = dataset("+2%",NFIN=0,range=(0,8),filename=filename3,directory='JCDijetBaseTask/jcdijet',color=2,style=24,rebin=4)
  
  JetScale_up = dataset("-2%",NFIN=0,range=(0,8),filename=filename4,directory='JCDijetBaseTask/jcdijet',color=2,style=24,rebin=4)
  #JetScale_down = dataset("-2%",NFIN=0,range=(0,8),filename=filename5,directory='JCDijetBaseTask/jcdijet',color=2,style=24,rebin=4)

  LHC13b = [Default_data,HadCorr_data]
  LHC13b_signal = [x.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = False) for x in LHC13b]
  jetPt = Default_data.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = True)[1]
  LHC13b_graphs = [x.getGraphs() for x in LHC13b]
  LHC13b_gausRMS = [x[0] for x in LHC13b_graphs]
  LHC13b_gammaRMS = [x[1] for x  in LHC13b_graphs]
  LHC13b_gausYield = [x[2] for x in LHC13b_graphs]
  LHC13b_gammaYield = [x[3] for x in LHC13b_graphs]
  LHC13b_systGausRMS_rel = defs.makeSystError(LHC13b_gausRMS[0],LHC13b_gausRMS[1],rel=True)
  LHC13b_systGausRMS_rel.SetName("GausRMS_hadcorrSystematics_relative")
  LHC13b_systGammaRMS_rel = defs.makeSystError(LHC13b_gammaRMS[0],LHC13b_gammaRMS[1],rel=True)
  LHC13b_systGammaRMS_rel.SetName("GammaRMS_hadcorrSystematics_relative")
  defs.drawErrors2(LHC13b_systGausRMS_rel,20,150,0,-0.10,0.10,0,"Narrow RMS","HadCorr","_data","SystematicErrors/SystematicErrorsGausRMS_EmcalHadCorr.pdf",0,0,0.02);
  defs.drawErrors2(LHC13b_systGammaRMS_rel,20,150,0,-0.10,0.10,0,"Wide RMS","HadCorr","_data","SystematicErrors/SystematicErrorsGammaRMS_EmcalHadCorr.pdf",0,0,0.02);

  colors = [1,2,3,4]
  names = ["Default","HadCorr"]

  n_figs = 8
  n_rows = 2
  fig, axs = defs.makegrid(4,n_figs//4,xlog=True,ylog=True,d=d,shareY=False,figsize=(15,7.5))
  axs = axs.reshape(8)
  axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone' + '\n Inclusive jT' + '\n',fontsize = 7)
  for signal,name,color,j in zip(LHC13b_signal,names,colors,range(10)):
    print("Plot {}".format(name))
    for jT,pT,ax,i in zip(signal[start:],jetPt[start:],axs[0:4],range(0,9)):
      print(i)
      jT.SetMarkerColor(color)
      jT.SetLineColor(color)
      plot = rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label=name,fmt='o',fillstyle='none',ecolor='blue') #Plot jT histogram, 
      ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
  
      ax.set_xlim([0.1,22]) #Set x-axis limits
      ax.set_ylim([5e-4,2e3]) #Set y-axis limits
      ax.set_xticklabels(ax.get_xticklabels(),horizontalalignment='left')
    ratios = []
    if j > 0:
      for jT,div in zip(signal,LHC13b_signal[0]):
        h = jT.Clone()
        h.Divide(div)
        ratios.append(h)
      axs[4].set_ylabel('Ratio',fontsize=9) #Add y-axis labels to left- and righmost subfigures
      axs[-1].set_ylabel('Ratio',fontsize=9)       
      for ratio,pT,ax in zip(ratios[start:],jetPt[start:],axs[4:9]):
        rplt.errorbar(ratio,xerr=False,emptybins=False,axes=ax,label='Ratio',fmt='o') #Plot ratio histogram,
        print("Draw {}".format(ratio.GetName()))   
        #if(i == 0):
        ax.set_yscale('linear')
        ax.set_xlim([0.1,20]) #Set x-axis limits
        ax.set_ylim([0.95,1.05]) #Set y-axis limits
     
    axs[0].legend(loc = 'lower left')
        
  plt.savefig("SystematicErrors/LHC13bHadCorrComparisonJetPt{}To{}.pdf".format(start,start+4),format='pdf') #Save figure
  plt.show() #Draw figure on screen
  
  return
      
  start = 4
  
  datasets = [Default,HadCorr_high,HadCorr_low,JetScale_up]
  signal,jetPt = Default.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = True)
  signal2 = HadCorr_high.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = False)
  signal3 = HadCorr_low.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = False)
  signal4 = JetScale_up.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = False)
  #signal5 = JetScale_down.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = False)

  signals = [signal,signal2,signal3,signal4]
  names = ["Default","-2%","+2%","Jet pT +2%"]
  if(True):
    graphs = [data.getGraphs() for data in datasets]
    gausRMS = [x[0] for x in graphs]
    gammaRMS = [x[1] for x  in graphs]
    gausYield = [x[2] for x in graphs]
    gammaYield = [x[3] for x in graphs]
    systGausRMS = defs.makeSystError(gausRMS[0],gausRMS[1])
    systGausRMS.SetName("GausRMS_emcalSystematics")
    systGammaRMS = defs.makeSystError(gammaRMS[0],gammaRMS[1])
    systGammaRMS.SetName("GammaRMS_emcalSystematics")
    systGausRMS_abs = defs.makeSystError(gausRMS[0],gausRMS[1],abs=True)
    systGausRMS_abs.SetName("GausRMS_emcalSystematics_absolute")
    systGammaRMS_abs = defs.makeSystError(gammaRMS[0],gammaRMS[1],abs=True)
    systGammaRMS_abs.SetName("GammaRMS_emcalSystematics_absolute")
    systGausRMS_rel = defs.makeSystError(gausRMS[0],gausRMS[1],rel=True)
    systGausRMS_rel.SetName("GausRMS_emcalSystematics_relative")
    systGausRMS_rel2 = defs.makeSystError(gausRMS[0],gausRMS[2],rel=True)
    systGausRMS_rel2.SetName("GausRMS_emcalSystematics_relative2")
    systGammaRMS_rel = defs.makeSystError(gammaRMS[0],gammaRMS[1],rel=True)
    systGammaRMS_rel.SetName("GammaRMS_emcalSystematics_relative")
    systGammaRMS_rel2 = defs.makeSystError(gammaRMS[0],gammaRMS[2],rel=True)
    systGammaRMS_rel2.SetName("GammaRMS_emcalSystematics_relative2")
    
    systGausRMS_rel3 = defs.makeSystError(gausRMS[0],gausRMS[3],rel=True)
    systGausRMS_rel3.SetName("GausRMS_emcalSystematics_relative")
    #systGausRMS_rel4 = defs.makeSystError(gausRMS[0],gausRMS[4],rel=True)
    #systGausRMS_rel4.SetName("GausRMS_emcalSystematics_relative2")
    systGammaRMS_rel3 = defs.makeSystError(gammaRMS[0],gammaRMS[3],rel=True)
    systGammaRMS_rel3.SetName("GammaRMS_emcalSystematics_relative")
    #systGammaRMS_rel4 = defs.makeSystError(gammaRMS[0],gammaRMS[4],rel=True)
    #systGammaRMS_rel4.SetName("GammaRMS_emcalSystematics_relative2")    
    
    
    defs.drawErrors2(systGausRMS_rel,20,150,0,-0.04,0.04,0,"Narrow RMS","-2%","_data","SystematicErrors/SystematicErrorsGausRMS_Emcal.pdf",0,0,0.01,error2=systGausRMS_rel2,title3="+2%");
    defs.drawErrors2(systGammaRMS_rel,20,150,0,-0.04,0.04,0,"Wide RMS","-2%","_data","SystematicErrors/SystematicErrorsGammaRMS_Emcal.pdf",0,0,0.01,error2=systGammaRMS_rel2,title3="+2%");

    
    defs.drawErrors2(systGausRMS_rel3,20,150,0,-0.04,0.04,0,"Narrow RMS","jet pT -2%","_data","SystematicErrors/SystematicErrorsGausRMS_jetScale.pdf",0,0,0.01);
    defs.drawErrors2(systGammaRMS_rel3,20,150,0,-0.04,0.04,0,"Wide RMS","jet pT -2%","_data","SystematicErrors/SystematicErrorsGammaRMS_jetScale.pdf",0,0,0.01);

        
    
  colors = [1,2,3,4]
  n_figs = 8
  n_rows = 2
  fig, axs = defs.makegrid(4,n_figs//4,xlog=True,ylog=True,d=d,shareY=False,figsize=(15,7.5))
  axs = axs.reshape(8)
  axs[1].text(0.12,0.002,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone' + '\n Inclusive jT' + '\n',fontsize = 7)
  for signal,name,color,j in zip(signals,names,colors,range(10)):
    print("Plot {}".format(name))
    for jT,pT,ax,i in zip(signal[start:],jetPt[start:],axs[0:4],range(0,9)):
      print(i)
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
        print("Draw {}".format(ratio.GetName()))   
        #if(i == 0):
        ax.set_yscale('linear')
        ax.set_xlim([0.1,20]) #Set x-axis limits
        ax.set_ylim([0.95,1.05]) #Set y-axis limits
     
    axs[0].legend(loc = 'lower left')
        
  plt.savefig("SystematicErrors/HadCorrComparisonJetPt{}To{}.pdf".format(start,start+4),format='pdf') #Save figure
  plt.show() #Draw figure on screen



  
if __name__ == "__main__": main()