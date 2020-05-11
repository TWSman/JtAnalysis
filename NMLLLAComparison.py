import logging
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
import rootpy
import defs
import re
import matplotlib
import ROOT
from ROOT import TF1,TMath,TGraphErrors
from rootpy.plotting import Canvas
from rootpy.plotting import Legend,Graph
from rootpy.io import root_open
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
from dataset import *
import math

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

def main(): 

  logjt_nmlla = [-0.8375,-0.68327,-0.50528,-0.32737,-0.16131,-0.042676,0.09963,0.23006,0.39612,0.56213,0.6451,0.81108,0.89409,1.0956,1.2022,1.3089,1.4274,1.5103,1.6169,1.7117,1.7945,1.8774,1.9601,2.0192,2.0901,2.1372,2.1723,2.1957,2.219,2.2187]
  nmlla = [5.6912,4.4526,3.4842,2.4101,1.7369,1.4155,1.0199,0.73479,0.52954,0.35149,0.27483,0.17508,0.14264,0.07711,0.053306,0.03685,0.025477,0.017609,0.011212,0.0074379,0.004735,0.0031408,0.0017674,0.0010796,0.00060746,0.00031475,0.0001502,0.000081077,0.000043764,0.00002672]
  jt_nmlla = [math.exp(logjt) for logjt in logjt_nmlla]
  nmlla = [y/x for x,y in zip(jt_nmlla,nmlla)]
  print(jt_nmlla)
  
  print('Number of arguments: ', len(sys.argv), 'arguments.')
  print('Argument list:',str(sys.argv))
  filename = sys.argv[1]
  filename2 = sys.argv[2]
  print("Input file data: ")
  print(filename)
  print("Input NMLLA: ")
  print(filename2)
  #separate = 5 #bin 5: 60-80 Gev
  Mixed_FullJets_R04 = datasetMixed("Full jets R=0.4",NFIN=0,range=5,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager',color=2,style=24,rebin=2)
  signal,jetPt = Mixed_FullJets_R04.getSubtracted('JetConeJtWeightBin','BgJtWeightBin',jetpt = True)
  for separate in range(2,8):
    Qs = {5:26,4:18,6:35,7:46,3:13,2:9}
    Q = Qs[separate] # 9,    13,   18,   26,   36,   46
    #Jet Pt: 22.5, 32.5, 45.0, 65.0, 87.5, 115.0]
    f_nmlla = root_open(filename2,"open")
    nmlla_fromfile = f_nmlla.get('NMLLAQ{:02d}'.format(Q))
    nmlla_gluon = f_nmlla.get('NMLLA_gluonQ{:02d}'.format(Q))
    nmlla_quark = f_nmlla.get('NMLLA_quarkQ{:02d}'.format(Q))
    scale = 1/nmlla_fromfile.Eval(1)
    print("Scale by {}".format(scale))
    nmlla_fromfile = scaleGraph(nmlla_fromfile,scale)
    nmlla_gluon = scaleGraph(nmlla_gluon,scale)
    nmlla_quark = scaleGraph(nmlla_quark,scale)
    nmlla_fromfile.SetMarkerColor(3)
    nmlla_quark.SetMarkerColor(2)
    nmlla_gluon.SetMarkerColor(3)
    nmlla_quark.SetLineColor(2)
    nmlla_gluon.SetLineColor(3)
    fit1,d1 = fitJtHisto(signal[separate],1,separate)
    fit2,d2 = fitJtHisto(nmlla_fromfile,1,separate)

    print(r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(jetPt[separate][0],jetPt[separate][1]))
    print('Data Wide RMS: {}, NMLLA Wide RMS: {}'.format(d1['gammaRMS'],d2['gammaRMS']))
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.text(0.2,0.0005,d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n Jet Cone',fontsize = 7)
    rplt.errorbar(signal[separate],xerr=False,emptybins=False,axes=ax,label=Mixed_FullJets_R04.name(),fmt='o') #Plot jT histogram,
    rplt.errorbar(nmlla_fromfile,xerr=False,yerr=False,emptybins=False,axes=ax,label="NMLLA",fmt='r-',markersize=0)
    #rplt.errorbar(nmlla_gluon,xerr=False,yerr=False,emptybins=False,axes=ax,label="NMLLA Gluon",fmt='r-')
    #rplt.errorbar(nmlla_quark,xerr=False,yerr=False,emptybins=False,axes=ax,label="NMLLA Quark",fmt='b-')
    if(separate==5):
      plt.plot(jt_nmlla,nmlla,axes=ax,label="NMLLA")
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(jetPt[separate][0],jetPt[separate][1])) 
    ax.set_xlim([0.1,12])
    ax.set_ylim([5e-6,1.5e3])
    ax.legend(loc = 'lower left')
    
    plt.savefig("PythonFigures/MixedFullJetsR04JetConeJtSignalJetPt{0}WithNMLLA.pdf".format(separate),format='pdf') #Save figure
    plt.show() #Draw figure on screen


def fitJtHisto(histo,cut,iJet):
  d = {}
  gaussfit = TF1("gaussfit","gausn",0,10)
  gaussfit.FixParameter(1,0)
  gaussfit.SetParLimits(2,B1low,B1high)
  gaussfit.SetParLimits(0,B2low,B2high)
  histo.Fit("gaussfit","QN")
  invG = TF1("invG","[0]* (pow([1],[2])/TMath::Gamma([2]))*exp(-[1]/x) * pow(x, -[2]-1)",0,10)
  invG.SetParameter(0,B3start[iJet])
  if iJet < B3cut:
    invG.SetParLimits(0,B3low,B3high1)
  else:
    invG.SetParLimits(0,B3low,B3high2)
  invG.SetParameter(1,B5start[iJet])
  invG.SetParLimits(1,B5low,B5high)
  
  invG.SetParameter(2,B4start[iJet])
  invG.SetParLimits(2,B4low,B4high)
  histo.Fit("invG","QN","",cut,3)
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
  gaussfit3.SetParLimits(4,B5low,B5high) #B5
  gaussfit3.SetParLimits(5,B4low,B4high) #B4
  try:
    lastBin = histo.FindLastBinAbove(1e-5)
    end = histo.GetBinCenter(lastBin)
  except AttributeError:
    end = 5
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
  return gaussfit3,d

def scaleGraph(gr1,sc):
  x1 = ROOT.Double()
  x_ = ROOT.Double()
  y1 = ROOT.Double()
  test = []
  y = []

  NC =  gr1.GetN()
  x = [0 for i in range(NC)]

  for ii in range(NC):
    gr1.GetPoint(ii, x_, y1)
    x1 = x_ * 1.0
    x[ii] = x1
    y.append(y1*sc)
    
  gr= Graph(NC)
  for x0,y0,i in zip(x,y,range(NC)):
    gr.SetPoint(i,x0,y0)
  return gr

  
if __name__ == "__main__": main()
