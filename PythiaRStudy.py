import logging

mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)
import rootpy
import defs
import ROOT
import re
import matplotlib
from rootpy.plotting import Canvas, Graph, Legend
from rootpy.io import root_open
from ROOT import TF1, TMath, TGraphErrors
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
from dataset import *
from array import array
from ctypes import c_int
import math


B1low = 0.1
B1start = [0.15, 0.15, 0.17, 0.18, 0.15, 0.15, 0.21, 0.15, 0.15]
B1high = 0.5
B2low = 40
B2high = 150
B2start = [50, 50, 66.20, 67.90, 50, 50, 75.89, 77.62, 77.62]
B3cut = 6
B3low = 0
B3high1 = 10
B3high2 = 25
B3start = [7, 7, 5.06, 4.90, 10, 10, 7, 12.93, 12.93]
B4low = 2.0
B4high = 15
B4start = [3, 3, 9.91, 8.80, 3, 3, 5.62, 4.18, 4.18]
B5low = 0.90
B5high = 4.5
B5start = [1.5, 1.5, 4.5, 4.5, 1.5, 1.5, 2.88, 1.62, 1.62]
peakstart = 0.5
peaklow = 0.3
peakhigh = 0.8
Rs = (0.3, 0.4, 0.5, -0.4)


def main():
    print("Number of arguments: ", len(sys.argv), "arguments.")
    print("Argument list:", str(sys.argv))
    filename = sys.argv[1]
    if len(sys.argv) > 2:
        fileData = sys.argv[2]
    else:
        fileData = None
    print("Input file: ")
    print(filename)
    Njets = 9
    # FullJets_R04 = dataset("FullR04",NFIN=6,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',color=2,style=24,rebin=5)
    # FullJets_R05 = dataset("FullR05",NFIN=7,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',color=3,style=24,rebin=5)
    # FullJets_R03 = dataset("FullR06",NFIN=8,filename=filename,directory='AliJJetJtTask/AliJJetJtHistManager',color=4,style=24,rebin=5)

    if fileData is not None:
        fD = root_open(fileData, "read")
        iS = 0
        gGausRMSData = fD.Get("gGausRMS{:02d}".format(iS))
        gGausRMSerrData = fD.Get("gGausRMS{:02d}_Systematics".format(iS))
        gGausYieldData = fD.Get("gGausYield{:02d}".format(iS))
        gGausYielderrData = fD.Get("gGausYield{:02d}_Systematics".format(iS))
        gGammaRMSData = fD.Get("gGammaRMS{:02d}".format(iS))
        gGammaRMSerrData = fD.Get("gGammaRMS{:02d}_Systematics".format(iS))
        gGammaYieldData = fD.Get("gGammaYield{:02d}".format(iS))
        gGammaYielderrData = fD.Get("gGammaYield{:02d}_Systematics".format(iS))

    # errGraph = [fD.Get('JetConeJtWeightBinNFin{:02d}JetPt{:02d}_Systematics'.format(iS,i)) for i in range(2,Njets)]  #Get jT histograms from file an array
    # hJtSignalGraph = [fD.Get('JetConeJtWeightBinNFin{:02d}JetPt{:02d}_Statistics'.format(iS,i)) for i in range(2,Njets)]  #Get jT histograms from file an array

    with root_open(filename, "read") as fin:
        FullJets_jT = []
        FullJets_fit = []
        FullJets_parameters = []
        FullJets_gausRMS = []
        FullJets_gammaRMS = []
        FullJets_gausYield = []
        FullJets_gammaYield = []
        colors = (1, 2, 3)
        for iF in (6, 7, 8):
            # iF = 8
            c = 1
            jT = [
                fin.get(
                    "AliJJetJtTask/AliJJetJtHistManager/JetConeJtWeightBin/JetConeJtWeightBinNFin{:02d}JetPt{:02d}".format(
                        iF, ij
                    )
                )
                for ij in range(Njets)
            ]
            bgJt = [
                fin.get(
                    "AliJJetJtTask/AliJJetJtHistManager/BgJtWeightBin/BgJtWeightBinNFin{:02d}JetPt{:02d}".format(
                        iF, ij
                    )
                )
                for ij in range(Njets)
            ]
            jetPtBin = [
                fin.get(
                    "AliJJetJtTask/AliJJetJtHistManager/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}".format(
                        iF, ij
                    )
                )
                for ij in range(Njets)
            ]
            nJets = [h.Integral() for h in jetPtBin]
            nBgs = [
                fin.get(
                    "AliJJetJtTask/AliJJetJtHistManager/BgTrkNumberBin/BgTrkNumberBinNFin{:02d}JetPt{:02d}".format(
                        iF, ij
                    )
                ).Integral()
                for ij in range(Njets)
            ]
            FullJets_jT.append(jT)
            for j, b, nj, nb in zip(jT, bgJt, nJets, nBgs):
                j.Rebin(4)
                b.Rebin(4)
                j.Scale(1.0 / nj, "width")
                b.Scale(1.0 / nb, "width")
                j.SetMarkerColor(c)
                b.SetMarkerColor(c)
                j.Add(b, -1.0)

        print(h.GetTitle())
        jetPt = [
            (
                int(
                    re.search(
                        r"p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*",
                        h.GetTitle(),
                        re.M | re.I,
                    ).group(1)
                ),
                int(
                    re.search(
                        r"p_{T,jet} : ([\d]*)\.[\d] - ([\d]*).[\d]*",
                        h.GetTitle(),
                        re.M | re.I,
                    ).group(2)
                ),
            )
            for h in FullJets_jT[0]
        ]  # Use regular expressions to extract jet pT range from histogram titles
        jetPtCenter = array("d", [(a + b) / 2.0 for a, b in jetPt])
        jetPtErrors = array("d", [(b - a) / 2.0 for a, b in jetPt])
        print(jetPt, type(jetPt))
        print(jetPtCenter, type(jetPtCenter))
        print(jetPtErrors, type(jetPtErrors))

        for i, jT in enumerate(FullJets_jT):
            print(i)
            print(jT)
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
            for h, i in zip(jT, range(Njets)):
                fit, d = fitJtHisto(h, "", 1, i, 8)
                fits.append(fit)
                parameters.append(d)
                gausRMS.append(d["gausRMS"])
                gausRMSe.append(d["gausRMSe"])
                gammaRMS.append(d["gammaRMS"])
                gammaRMSe.append(d["gammaRMSe"])
                gausYield.append(d["gausYield"])
                gausYielde.append(d["gausYielde"])
                gammaYield.append(d["gammaYield"])
                gammaYielde.append(d["gammaYielde"])

            gausRMSg = Graph(len(gausRMS) - 2)
            gammaRMSg = Graph(len(gammaRMS) - 2)
            gausYieldg = Graph(len(gausYield) - 2)
            gammaYieldg = Graph(len(gammaYield) - 2)
            for h, he, g in zip(
                (gausYield, gammaYield),
                (gausYielde, gammaYielde),
                (gausYieldg, gammaYieldg),
            ):
                for x, xe, a, e, i in zip(
                    jetPtCenter[2:],
                    jetPtErrors[2:],
                    h[2:],
                    he[2:],
                    range(len(gausRMS) - 2),
                ):
                    g.SetPoint(i, x, a)
                    g.SetPointError(i, xe, xe, e, e)

            for a, b, c, d, e, f, i in zip(
                gausRMS[2:],
                gammaRMS[2:],
                gausRMSe[2:],
                gammaRMSe[2:],
                jetPtCenter[2:],
                jetPtErrors[2:],
                range(len(gausRMS) - 2),
            ):
                gausRMSg.SetPoint(i, e, a)
                gausRMSg.SetPointError(i, f, f, c, c)
                gammaRMSg.SetPoint(i, e, b)
                gammaRMSg.SetPointError(i, f, f, d, d)

            FullJets_gausRMS.append(gausRMSg)
            print(len(FullJets_gausRMS))
            FullJets_gammaRMS.append(gammaRMSg)
            FullJets_gausYield.append(gausYieldg)
            FullJets_gammaYield.append(gammaYieldg)
            FullJets_fit.append(fits)
            FullJets_parameters.append(parameters)

        if fileData is not None:

            FullJets_gausRMS.append(gGausRMSData)
            FullJets_gammaRMS.append(gGammaRMSData)
            FullJets_gausYield.append(gGausYieldData)
            FullJets_gammaYield.append(gGammaYieldData)

        print(len(FullJets_gausRMS))
        drawWithErrors2Combined(
            FullJets_gausRMS,
            FullJets_gammaRMS,
            15,
            500,
            1,
            0,
            1.65,
            0,
            r"$p_\mathrm{T,jet}$",
            r"$\sqrt{\left<j_\mathrm{T}^2\right>}$",
            "Pythia Truth",
            "PythonFigures/RcomparisonRMS",
        )
        drawWithErrors2Combined(
            FullJets_gausYield,
            FullJets_gammaYield,
            15,
            500,
            1,
            0,
            10,
            0,
            r"$p_\mathrm{T,jet}$",
            r"Yield",
            "Pythia Truth",
            "PythonFigures/RcomparisonYield",
        )

        ratios = []
        for hists in FullJets_jT:
            if hists is not None:
                ratio = []
                for h, true in zip(hists, FullJets_jT[1]):
                    h2 = h.Clone()
                    h2.Divide(true)
                    ratio.append(h2)
            else:
                ratio = None
            ratios.append(ratio)

        fig, axs = defs.makegrid(4, 2, xlog=True, ylog=False, d=d, shareY=False)
        axs = axs.reshape(8)
        axs[4].set_ylabel("Ratio to R = 0.4")
        axs[7].set_ylabel("Ratio to R = 0.4")

        axs[1].text(
            0.02,
            0.005,
            "Pythia\n"
            r"pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$"
            "\n Full jets \n"
            r"Anti-$k_\mathrm{T}$",
            fontsize=7,
        )  # Add text to second subfigure, first parameters are coordinates in the drawn scale/units
        for hists, R in zip(FullJets_jT, Rs):
            for jT, ax, i, pT in zip(hists[2:], axs[0:4], range(0, 9), jetPt[2:]):
                rplt.errorbar(
                    jT,
                    emptybins=False,
                    xerr=False,
                    label="R = {:.1f}".format(R),
                    axes=ax,
                    fmt="o",
                )  # Plot jT histogram,
                ax.text(
                    0.3,
                    1e2,
                    r"$p_\mathrm{T,jet}$:"
                    "\n"
                    r" {:02d}-{:02d} GeV".format(pT[0], pT[1]),
                )
                ax.set_xlim([0.01, 20])  # Set x-axis limits
                ax.set_ylim([5e-4, 2e3])  # Set y-axis limits
                ax.set_yscale("log")
                ax.grid(True)

        for ratio in ratios:
            for r, ax in zip(ratio[2:], axs[4:8]):
                rplt.errorbar(r, axes=ax, emptybins=False, xerr=False)
                ax.set_xlim([0.01, 20])
                ax.set_ylim([0.1, 3])
                ax.grid(True)

        axs[0].legend(loc="lower left")

        plt.savefig("PythonFigures/RcomparisonSignal.pdf", format="pdf")  # Save figure
        plt.show()  # Draw figure on screen


def drawWithErrors2Combined(
    hists, hists2, xlow, xhigh, logx, ylow, yhigh, ylog, xTitle, yTitle, title, file
):
    fig, axs = plt.subplots(2, 2, figsize=(7, 7))
    axs = axs.reshape(4)
    ax = axs[0]
    ax2 = axs[1]
    ax.set_xlim([xlow, xhigh])
    ax.set_ylim([ylow, yhigh])
    ax.set_xlabel(xTitle)  # Add x-axis labels for bottom row
    ax.set_ylabel(yTitle)  # Add x-axis labels for bottom row
    ratios1 = []
    ratios2 = []
    if logx:
        ax.set_xscale("log")
        ax2.set_xscale("log")
    for h, h2, c, R in zip(hists, hists2, (1, 2, 3, 4), Rs):
        h.SetMarkerColor(c)
        h.SetMarkerStyle(24)
        h2.SetMarkerColor(c)
        h2.SetMarkerStyle(25)
        if len(hists) > 1:
            ratios1.append(grrDivide(h, hists[1]))
            ratios2.append(grrDivide(h2, hists2[1]))
        if c == 1:
            rplt.errorbar(
                h, xerr=False, emptybins=False, label="Narrow", axes=ax, fmt="+"
            )
        else:
            rplt.errorbar(h, xerr=False, emptybins=False, axes=ax, fmt="+")
        # rplt.errorbar(h,xerr=False,emptybins=False,axes=ax,fmt='+')

        if R > 0:
            rplt.errorbar(
                h2,
                xerr=False,
                emptybins=False,
                axes=ax2,
                label="R = {:.1f} Wide".format(R),
                fmt="+",
            )
        else:
            rplt.errorbar(
                h2,
                xerr=False,
                emptybins=False,
                axes=ax2,
                label="R = {:.1f} Data".format(-R),
                fmt="+",
            )

    print(
        "({} - {})/9.0 + {} is {}".format(
            yhigh, ylow, ylow, (yhigh - ylow) / 9.0 + ylow
        )
    )
    ax.text(
        100,
        (yhigh - ylow) / 3.0 + ylow,
        title
        + "\n"
        + d["system"]
        + "\n"
        + d["jettype"]
        + "\n"
        + r"Anti-$k_\mathrm{T}$",
        fontsize=10,
    )

    ax.legend(loc="upper left")
    ax.set_xlim([xlow, xhigh])
    ax.set_ylim([ylow, yhigh])
    ax.grid(True)
    ax2.set_xlim([xlow, xhigh])
    ax2.set_ylim([ylow, yhigh])
    ax2.grid(True)

    ax = axs[2]
    ax.grid(True)
    ax.set_xlabel(xTitle)  # Add x-axis labels for bottom row
    ax.set_ylabel("Ratio")  # Add x-axis labels for bottom row
    ax2 = axs[3]
    ax2.grid(True)
    ax2.set_xlabel(xTitle)  # Add x-axis labels for bottom row
    ax2.set_ylabel("Ratio")  # Add x-axis labels for bottom row
    for ratio, ratio2, c in zip(ratios1, ratios2, (1, 2, 3, 4)):
        ratio.Shift(-2)
        ratio2.Shift(2)
        ratio.SetMarkerColor(c)
        ratio2.SetMarkerColor(c)
        ratio.SetMarkerStyle(24)
        ratio2.SetMarkerStyle(25)
        rplt.errorbar(ratio, xerr=False, emptybins=False, axes=ax)
        rplt.errorbar(ratio2, xerr=False, emptybins=False, axes=ax2)
    ax.set_xlim([xlow, xhigh])
    ax.set_ylim([0.8, 1.22])
    ax2.set_xlim([xlow, xhigh])
    ax2.set_ylim([0.8, 1.22])
    if logx:
        ax.set_xscale("log")
        ax2.set_xscale("log")

    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)  # Set space between subfigures to 0
    plt.savefig("{}.pdf".format(file), format="pdf")  # Save figure

    plt.show()  # Draw figure on screen


def grrDivide(gr1, gr2):
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

    NC = gr1.GetN()
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
        y.append(y1 / y2 if y2 != 0 else 0)
        ex.append(x2e)
        ey.append(
            math.sqrt(pow(y1e / y2, 2) + pow(y1 * y2e / (y2 * y2), 2)) if y2 != 0 else 0
        )

    gr = Graph(NC)
    for x0, y0, x0e, y0e, i in zip(x, y, ex, ey, range(NC)):
        gr.SetPoint(i, x0, y0)
        gr.SetPointError(i, x0e, x0e, y0e, y0e)
    return gr


def fitJtHisto(histo, method, cut, iJet, iFinder):
    d = {}
    gaussfit = TF1("gaussfit", "gausn", 0, 10)
    gaussfit.FixParameter(1, 0)
    gaussfit.SetParLimits(2, B1low, B1high)
    gaussfit.SetParLimits(0, B2low, B2high)
    histo.Fit("gaussfit", "QN")
    if method == "alt":
        invG = TF1(
            "invG",
            "[0]* (pow([1]*([2]+1),[2])/TMath::Gamma([2]))*exp(-[1]*([2]+1)/x) * pow(x, -[2]-1)",
            0,
            10,
        )
    else:
        invG = TF1(
            "invG",
            "[0]* (pow([1],[2])/TMath::Gamma([2]))*exp(-[1]/x) * pow(x, -[2]-1)",
            0,
            10,
        )
    invG.SetParameter(0, B3start[iJet])
    if iJet < B3cut:
        invG.SetParLimits(0, B3low, B3high1)
    else:
        invG.SetParLimits(0, B3low, B3high2)
    if method == "alt":
        invG.SetParameter(1, peakstart)
        invG.SetParLimits(1, peaklow, peakhigh)
    else:
        invG.SetParameter(1, B5start[iJet])
        invG.SetParLimits(1, B5low, B5high)

    invG.SetParameter(2, B4start[iJet])
    invG.SetParLimits(2, B4low, B4high)
    histo.Fit("invG", "QN", "", cut, 3)
    if method == "alt":
        gaussfit3 = TF1(
            "gaussfit3",
            "gausn(0) + [3]* (pow([4]*([5]+1),[5])/TMath::Gamma([5]))*exp(-[4]*([5]+1)/x) * pow(x, -[5]-1)",
            0,
            10,
        )
    else:
        gaussfit3 = TF1(
            "gaussfit3",
            "gausn(0) + [3]* (pow([4],[5])/TMath::Gamma([5]))*exp(-[4]/x) * pow(x, -[5]-1)",
            0,
            10,
        )
    gaussfit3.SetParameter(0, gaussfit.GetParameter(0))
    gaussfit3.SetParameter(1, gaussfit.GetParameter(1))
    gaussfit3.SetParameter(2, gaussfit.GetParameter(2))
    gaussfit3.SetParameter(3, invG.GetParameter(0))
    gaussfit3.SetParameter(4, invG.GetParameter(1))
    gaussfit3.SetParameter(5, invG.GetParameter(2))
    gaussfit3.FixParameter(1, 0)

    gaussfit3.SetParLimits(0, B2low, B2high)  # B2
    gaussfit3.SetParLimits(2, B1low, B1high)  # B1
    if iJet < B3cut:
        gaussfit3.SetParLimits(3, B3low, B3high1)  # B3
    else:
        gaussfit3.SetParLimits(3, B3low, B3high2)  # B3
    if method == "alt":
        gaussfit3.SetParLimits(4, peaklow, peakhigh)  # Peak
    else:
        gaussfit3.SetParLimits(4, B5low, B5high)  # B5
    gaussfit3.SetParLimits(5, B4low, B4high)  # B4

    lastBin = histo.FindLastBinAbove(1e-5)
    end = histo.GetBinCenter(lastBin)
    histo.Fit("gaussfit3", "QN", "", 0.1, end)
    # chi2[iF][ij] = gaussfit3.GetChisquare()
    # chi2[iF][ij] = getChi2(histo,gaussfit3,0.1,end,0)
    # chi2dof[iF][ij] = chi2[iF][ij]/(bins[iF][ij]-5)
    # chi2n[iF][ij] = chi2[iF][ij]/bins[iF][ij]
    B2 = gaussfit3.GetParameter(0)
    B1 = gaussfit3.GetParameter(2)
    B2e = gaussfit3.GetParError(0)
    B1e = gaussfit3.GetParError(2)
    gaussigma = TMath.Sqrt(2) * B1
    gausyield = B2 * B1 / TMath.Sqrt(2 * TMath.Pi())
    gausyielde = TMath.Sqrt(
        (B2 * B2 * B1e * B1e + B1 * B1 * B2e * B2e) / (2 * TMath.Pi())
    )
    gaussigmae = TMath.Sqrt(2) * B1e

    constant = gaussfit3.GetParameter(3)  # B3
    constante = gaussfit3.GetParError(3)
    alpha = gaussfit3.GetParameter(5)  # B4
    alphae = gaussfit3.GetParError(5)
    if method == "alt":
        peak = gaussfit3.GetParameter(4)
        peake = gaussfit3.GetParError(4)
        beta = peak * (alpha + 1)  # B5
        betae = TMath.Sqrt(pow((alpha + 1) * peake, 2) + pow(peak * alphae, 2))
    else:
        beta = gaussfit3.GetParameter(4)  # B5
        betae = gaussfit3.GetParError(4)
        peak = beta / (alpha + 1)
        peake = TMath.Sqrt(
            pow(betae / (alpha + 1), 2) + pow(alphae * beta / pow(alpha + 1, 2), 2)
        )

    gammaYield = constant * beta / (alpha - 1)
    gammaRMS = beta / TMath.Sqrt((alpha - 2) * (alpha - 3))
    gammaYielde = TMath.Sqrt(
        pow(beta * constante / (alpha - 1), 2)
        + pow(constant * beta * alphae / pow(alpha - 1, 2), 2)
        + pow(constant * betae / (alpha - 1), 2)
    )
    gammaRMSe = TMath.Sqrt(
        pow(
            (5 - 2 * alpha) * beta * alphae / pow(2 * ((alpha - 2) * (alpha - 3)), 1.5),
            2,
        )
        + pow(betae / TMath.Sqrt((alpha - 2) * (alpha - 3)), 2)
    )

    d["B2"] = B2
    d["B1"] = B1
    d["B3"] = constant
    d["B4"] = alpha
    d["B5"] = beta
    d["peak"] = peak
    d["B2e"] = B2e
    d["B1e"] = B1e
    d["B3e"] = constante
    d["B4e"] = alphae
    d["B5e"] = betae
    d["peake"] = peake

    d["gausRMS"] = gaussigma
    d["gausRMSe"] = gaussigmae
    d["gausYield"] = gausyield
    d["gausYielde"] = gausyielde
    d["gammaRMS"] = gammaRMS
    d["gammaRMSe"] = gammaRMSe
    d["gammaYield"] = gammaYield
    d["gammaYielde"] = gammaYielde
    return gaussfit3, d


if __name__ == "__main__":
    main()
