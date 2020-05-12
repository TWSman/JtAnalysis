#!/usr/bin/python3
# dataset.py by Tomas Snellman
# Class to hold one dataset in Jet jT analysis


import rootpy
import defs
import re
from rootpy.plotting import Graph
from rootpy.io import root_open
from ROOT import TF1
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import sys
import math


def parse_jet_pt_bins(hist_list, search="jet"):
    if search == "jet":
        search_string = r"p_{T, ?jet} : ([\d]*)\.?[\d]? - ([\d]*)\.?[\d]*"
    else:
        search_string = r"p_{T, ?constituent}:([\d]*)\.?[\d]?-([\d]*)\.?[\d]*"

    # Use regular expressions to extract jet pT range from histogram titles
    res = [re.search(search_string, h.GetTitle(), re.M | re.I,) for h in hist_list]

    bins = [(int(re.group(1)), int(re.group(2))) for re in res]
    print(bins)
    return bins


class dataset(object):
    """
    Basic class for handling a single dataset with one input file/directory and one jet finder

    Can retrieve histograms binned by jet Pt.

    Author Tomas Snellman tomas.snellman@cern.ch

    Attributes:
      properties: List of kwargs used to create class
      _filename: Name of used file
      _NFIN: Index of jet finder
      _name: name of class. Used in figure legends
      _directory: Directory in the file where histograms are located
    """

    def __init__(self, name, **kwargs):
        """
        The constructor,

        Args:
        name: The name of the dataset, will show in figure legends

        Kwargs:
        **filename** Name of the file\n
        **directory** Directory where the histograms are to be found\n
        **NFIN** Jet finder index\n
        **rebin** How much rebinning should be done\n
        **range** Which histograms should be retrieved, 9 by default which gives the full jet pT range from 5 -10GeV bin to 150-150GeV bin\n
        """
        print("Creating {}".format(name))
        self._NFIN = kwargs["NFIN"]
        self._filename = kwargs["filename"]
        self._name = name
        self._directory = kwargs["directory"]
        self.properties = kwargs
        self._f = root_open(self._filename)
        self._rebin = kwargs.get("rebin", 1)
        self._measN = None
        self._measBgN = None
        self._measRndmBgN = None
        self._range = kwargs.get("range", (0, 9))
        self._fitDone = False
        self._drawFit = False

        with root_open(self._filename, "read") as f:
            self._measN = [
                f.Get(
                    "{}/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}".format(
                        self._directory, self._NFIN, i
                    )
                ).Integral()
                for i in range(self._range[0], self._range[1])
            ]  # Get number of jets by jet pT bins
            self._measBgN = [
                f.Get(
                    "{}/BgTrkNumberBin/BgTrkNumberBinNFin{:02d}JetPt{:02d}".format(
                        self._directory, self._NFIN, i
                    )
                ).Integral()
                for i in range(self._range[0], self._range[1])
            ]  # Get number of background jets
            try:
                self._measRndmBgN = [
                    int(
                        f.Get(
                            "{}/hNumber/hNumberNFin{:02d}".format(
                                self._directory, self._NFIN
                            )
                        ).GetBinContent(7 + i)
                    )
                    for i in range(self._range[0], self._range[1])
                ]  # Get number of background jets
            except rootpy.io.file.DoesNotExist:
                print("No random background")
                self._measRndmBgN = [0 for i in range(self._range[0], self._range[1])]

    def printStats(self, **kwargs):
        """
        Print available statistics, number of jets by jet PT bin

        Args:

        Kwargs:
        **format** if set to latex will print a latex table source code, otherwise python lists \n
        **what** What statistics to give \n
        jets: gives number of jets\n
        bg: gives number of background cones \n
        bgratio: Gives ratio of bg cones to number of jets \n
        """
        ratios = [a / float(b) for b, a in zip(self._measN, self._measBgN)]
        hist = [
            self._f.Get(
                "{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}".format(
                    {
                        "dir": self._directory,
                        "histname": "JetPtBin",
                        "NFin": self._NFIN,
                        "pT": i,
                    }
                )
            )
            for i in range(self._range[0], self._range[1])
        ]  # Get jT histograms from file an array
        jetPt = parse_jet_pt_bins(hist)

    def getSubtracted(self, inclusive, background, **kwargs):
        hist, jetpt = self.getHist(inclusive, jetpt=True)
        if kwargs.get("randomBG", False):
            bg = self.getHist(background, isRndmBg=True, jetpt=False)
        else:
            bg = self.getHist(background, isBg=True, jetpt=False)
        signals = []
        for h, b in zip(hist, bg):
            s = h.Clone()
            s.Add(b, -1)
            signals.append(s)
        if kwargs.get("jetpt", False):
            return signals, jetpt
        else:
            return signals

    def getHistogram(self, name, default, **kwargs):
        h = self._f.get(name)
        if h is not None:
            return h
        else:
            h = default.Clone()
            h.clear()
            return h

    def getJetPt(self, **kwargs):
        try:
            if "dir" in kwargs:
                print(
                    "{0[dir]}/JetPt/JetPtNFin{0[NFin]:02d}".format(
                        {"dir": kwargs["dir"], "NFin": self._NFIN}
                    )
                )
                hist = self._f.Get(
                    "{0[dir]}/JetPt/JetPtNFin{0[NFin]:02d}".format(
                        {"dir": kwargs["dir"], "NFin": self._NFIN}
                    )
                )
            else:
                print(
                    "{0[dir]}/JetPt/JetPtNFin{0[NFin]:02d}".format(
                        {"dir": self._directory, "NFin": self._NFIN}
                    )
                )
                hist = self._f.Get(
                    "{0[dir]}/JetPt/JetPtNFin{0[NFin]:02d}".format(
                        {"dir": self._directory, "NFin": self._NFIN}
                    )
                )
        except rootpy.io.file.DoesNotExist:
            return None
        hist.Scale(1.0, "width")
        return hist

    def setDrawFit(self, drawF):
        self._drawFit = drawF

    def getHistByName(self, name, ij, **kwargs):
        """
          Retrieve a single histogram by name
        """
        try:
            if "dir" in kwargs:
                hist = self._f.Get(
                    "{0[dir]}/{0[histname]}".format(
                        {"dir": kwargs["dir"], "histname": name}
                    )
                ).Clone()  # Get a single jT histogram from file
            else:
                hist = self._f.Get(
                    "{0[dir]}/{0[histname]}".format(
                        {"dir": self._directory, "histname": name}
                    )
                ).Clone()  # Get a single jT histogram from file
        except rootpy.io.file.DoesNotExist:
            return None
        hist.Sumw2()
        hist.Rebin(self._rebin)
        if kwargs.get("isBg", False):
            N = self._measBgN[ij]
        else:
            N = self._measN[ij]
        if self.properties.get("isWeight", False) or kwargs.get("isWeight", False):
            hist.Scale(1.0, "width")
            print("ISWEIGHT")
        else:
            hist.Scale(1.0 / N, "width")
        return hist

    def getHist(self, name, **kwargs):
        """
        Retrieve a list of histograms by jet pT bins

        Args:
          name: Name of histogram

        Kwargs:
          isbg: Determines which normalization to use
          isWeight: Determines if the data is already normalized by number of jets/background. Bin width correction is done anyway
          extra: String to add to the end of the histogram name. (For additional binning etc.)

        """
        extra = kwargs.get("extra", "")
        print(name)
        print(extra)
        format_string = "{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}{0[extra]}"

        if "dir" in kwargs:
            hist = [
                self._f.Get(
                    format_string.format(
                        {
                            "dir": kwargs["dir"],
                            "histname": name,
                            "NFin": self._NFIN,
                            "pT": i,
                            "extra": extra,
                        }
                    )
                ).Clone()
                for i in range(self._range[0], self._range[1])
            ]  # Get jT histograms from file an array
        else:
            hist = [
                self._f.Get(
                    format_string.format(
                        {
                            "dir": self._directory,
                            "histname": name,
                            "NFin": self._NFIN,
                            "pT": i,
                            "extra": extra,
                        }
                    )
                ).Clone()
                for i in range(self._range[0], self._range[1])
            ]  # Get jT histograms from file an array
        # print('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory, 'histname':name,'NFin':self._NFIN,'pT':1}))
        jetPt = parse_jet_pt_bins(hist)

        if "LeadingRef" in name:
            normalization = [
                self._f.Get(
                    "{}/LeadingRefJetPtBin/LeadingRefJetPtBinNFin{:02d}JetPt{:02d}".format(
                        self._directory, self._NFIN, i
                    )
                ).Integral()
                for i in range(self._range[0], self._range[1])
            ]  # Get number of jets by jet pT bins
            print("Normalization set to LeadingRef")
            print(normalization)
            print("Before:")
            print(self._measN)
        else:
            normalization = self._measN
        for h, N, bgN, rndmbgN in zip(
            hist, normalization, self._measBgN, self._measRndmBgN
        ):
            h.Sumw2()
            # print("Rebinning {} by {} in set {} that has {} bins".format(h.GetTitle(), self._rebin, self._name, h.GetNbinsX()))
            h.Rebin(self._rebin)
            print(kwargs)
            if self.properties.get("isWeight", False):
                h.SetLineColor(self.properties.get("color", 1))
                h.SetMarkerColor(self.properties.get("color", 1))
                h.Scale(1.0, "width")
            else:
                if kwargs.get("isBg", False):
                    h.SetLineColor(self.properties.get("color", 1) + 1)
                    h.SetMarkerColor(self.properties.get("color", 1) + 1)
                    h.Scale(1.0 / bgN, "width")
                    print("{} is bg".format(name))
                elif kwargs.get("isRndmBg", False):
                    print("Is random background")
                    h.SetLineColor(self.properties.get("color", 1) + 2)
                    h.SetMarkerColor(self.properties.get("color", 1) + 2)
                    h.Scale(1.0 / rndmbgN, "width")
                else:
                    h.SetLineColor(self.properties.get("color", 1))
                    h.SetMarkerColor(self.properties.get("color", 1))
                    h.Scale(1.0 / N, "width")

            h.SetMarkerStyle(self.properties.get("style", 24))
            h.SetMarkerSize(0.5)
            h.SetLineColor(1)

        if kwargs.get("jetpt", False):
            return hist, jetPt
        else:
            return hist

    def getHistPta(self, name, **kwargs):
        """
        Retrieve a list of histograms by jet pT bins and track pT bins

        Args:
          name: Name of histogram

        Kwargs:
          isbg: Determines which normalization to use
          isWeight: Determines if the data is already normalized by number of jets/background. Bin width correction is done anyway

        """
        hists = []
        format_string = "{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[jetpT]:02d}TrkPt{0[trkpT]:02d}"

        for i in range(self._range[0], self._range[1]):
            if "dir" in kwargs:
                hist = [
                    self._f.Get(
                        format_string.format(
                            {
                                "dir": kwargs["dir"],
                                "histname": name,
                                "NFin": self._NFIN,
                                "jetpT": i,
                                "trkpT": j,
                            }
                        )
                    ).Clone()
                    for j in range(0, 11)
                ]  # Get jT histograms from file an array
            else:
                hist = [
                    self._f.Get(
                        format_string.format(
                            {
                                "dir": self._directory,
                                "histname": name,
                                "NFin": self._NFIN,
                                "jetpT": i,
                                "trkpT": j,
                            }
                        )
                    ).Clone()
                    for j in range(0, 11)
                ]  # Get jT histograms from file an array
            hists.append(hist)
        # print('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory, 'histname':name,'NFin':self._NFIN,'pT':1}))

        # Get Jet Pt bins
        jetPt = parse_jet_pt_bins(hist)

        # Get Track pt Bins
        trkPt = parse_jet_pt_bins(search="constituent")

        # print(len(hist))
        # print(hist)
        # print(jetPt)
        for hist, N, bgN in zip(hists, self._measN, self._measBgN):
            for h in hist:
                h.Sumw2()
                # print("Rebinning {} by {} in set {} that has {} bins".format(h.GetTitle(), self._rebin, self._name, h.GetNbinsX()))
                h.Rebin(self._rebin)
                print(kwargs)
                if self.properties.get("isWeight", False):
                    h.SetLineColor(self.properties.get("color", 1))
                    h.SetMarkerColor(self.properties.get("color", 1))
                    h.Scale(1.0, "width")
                else:
                    if kwargs.get("isBg", False):
                        h.SetLineColor(self.properties.get("color", 1) + 1)
                        h.SetMarkerColor(self.properties.get("color", 1) + 1)
                        h.Scale(1.0 / bgN, "width")
                        print("{} is bg".format(name))
                    else:
                        h.SetLineColor(self.properties.get("color", 1))
                        h.SetMarkerColor(self.properties.get("color", 1))
                        h.Scale(1.0 / N, "width")

                h.SetMarkerStyle(self.properties.get("style", 24))
                h.SetMarkerSize(0.5)
                h.SetLineColor(1)

        if kwargs.get("jetpt", False):
            return hist, jetPt, trkPt
        else:
            return hist

    def get2DHist(self, name, **kwargs):
        """
        Retrieve a list of 2D histograms by jet pT bins

        Args:
          name: Name of histogram

        Kwargs:
          isbg: Determines which normalization to use

        """
        format_string = "{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}"

        if "dir" in kwargs:
            hist = [
                self._f.Get(format_string.format(
                        {
                            "dir": kwargs["dir"],
                            "histname": name,
                            "NFin": self._NFIN,
                            "pT": i,
                        }
                    )
                ).Clone()
                for i in range(self._range[0], self._range[1])
            ]  # Get jT histograms from file an array
        else:
            hist = [
                self._f.Get(format_string.format(
                        {
                            "dir": self._directory,
                            "histname": name,
                            "NFin": self._NFIN,
                            "pT": i,
                        }
                    )
                ).Clone()
                for i in range(self._range[0], self._range[1])
            ]  # Get jT histograms from file an array
        jetPt = parse_jet_pt_bins(hist)

        if kwargs.get("jetpt", False):
            return hist, jetPt
        else:
            return hist

    def getGraphs(self, **kwargs):
        return self.getGraphs2(self._range[0], self._range[1], **kwargs)

    def getGraphs2(self, start, end, **kwargs):
        if not self._fitDone:
            print("Fit hasn't been performed. Fitting now.")
            self._fitDone = True
            self._gausRMS = []
            self._gammaRMS = []
            self._gausRMSe = []
            self._gammaRMSe = []
            self._gausYield = []
            self._gammaYield = []
            self._gausYielde = []
            self._gammaYielde = []
            self._fits = []
            self._parameters = []
            jt, jetpt = self.getSubtracted(
                "JetConeJtWeightBin", "BgJtWeightBin", jetpt=True
            )
            print(jetpt)
            jetPtCenter = [(a + b) / 2 for a, b in jetpt]
            jetPtErrors = [(a - b) / 2 for a, b in jetpt]
            print(jetPtCenter)
            for h, iJet in zip(jt, range(start, end)):
                fit, d = fitJtHisto(h, "", 1, iJet, self._NFIN, drawFit=self._drawFit)
                self._fits.append(fit)
                self._parameters.append(d)
                self._gausRMS.append(d["gausRMS"])
                self._gausRMSe.append(d["gausRMSe"])
                self._gammaRMS.append(d["gammaRMS"])
                self._gammaRMSe.append(d["gammaRMSe"])
                self._gausYield.append(d["gausYield"])
                self._gausYielde.append(d["gausYielde"])
                self._gammaYield.append(d["gammaYield"])
                self._gammaYielde.append(d["gammaYielde"])
            self._gausRMSg = Graph(len(self._gausRMS) - 2)
            self._gammaRMSg = Graph(len(self._gammaRMS) - 2)
            self._gausYieldg = Graph(len(self._gausYield) - 2)
            self._gammaYieldg = Graph(len(self._gammaYield) - 2)
            for h, he, g in zip(
                (self._gausYield, self._gammaYield),
                (self._gausYielde, self._gammaYielde),
                (self._gausYieldg, self._gammaYieldg),
            ):
                for x, xe, a, e, i in zip(
                    jetPtCenter[2:],
                    jetPtErrors[2:],
                    h[2:],
                    he[2:],
                    range(len(self._gausRMS) - 2),
                ):
                    g.SetPoint(i, x, a)
                    g.SetPointError(i, xe, xe, e, e)

            for a, b, c, d, e, f, i in zip(
                self._gausRMS[2:],
                self._gammaRMS[2:],
                self._gausRMSe[2:],
                self._gammaRMSe[2:],
                jetPtCenter[2:],
                jetPtErrors[2:],
                range(len(self._gausRMS) - 2),
            ):
                self._gausRMSg.SetPoint(i, e, a)
                self._gausRMSg.SetPointError(i, f, f, c, c)
                self._gammaRMSg.SetPoint(i, e, b)
                self._gammaRMSg.SetPointError(i, f, f, d, d)

            for a, b, c, d, e, f, i in zip(
                self._gausYield[2:],
                self._gammaYield[2:],
                self._gausYielde[2:],
                self._gammaYielde[2:],
                jetPtCenter[2:],
                jetPtErrors[2:],
                range(len(self._gausYield) - 2),
            ):
                self._gausYieldg.SetPoint(i, e, a)
                self._gausYieldg.SetPointError(i, f, f, c, c)
                self._gammaYieldg.SetPoint(i, e, b)
                self._gammaYieldg.SetPointError(i, f, f, d, d)

        return (self._gausRMSg, self._gammaRMSg, self._gausYieldg, self._gammaYieldg)

    def name(self):
        return self._name

    def getprop(self, key):
        return self.properties.get(key, None)


d = dict(
    jetalg=r"Anti-$k_\mathrm{T}$, R=0.4",
    jettype="Full Jets",
    system=r"ALICE p-Pb " + "\n" + r"$\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$",
    trigger="kINT7/kEMCEJE",
    cut=r"$\left| \eta_\mathrm{jet} \right| < 0.25$",
)


class datasetMixed(dataset, object):
    """Class for handling datasets where different bins come from
    different files.
    For example MB data for low jet pT and triggered data for high jet pT

    """

    def __init__(self, name, **kwargs):
        """Constructor

        """
        super(datasetMixed, self).__init__(name, **kwargs)
        self._filename1 = kwargs["filename"]
        self._filename2 = kwargs.get("filename2", self._filename1)
        if self._filename2 != self._filename1:
            self._f2 = root_open(self._filename)
        else:
            self._f2 = self._f
        self._directory1 = kwargs["directory"]
        self._directory2 = kwargs["directory2"]
        print(
            "Directory 1: {} Directory 2: {}".format(self._directory1, self._directory2)
        )
        self._measN1 = [
            self._f.Get(
                "{}/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}".format(
                    self._directory1, self._NFIN, i
                )
            ).Integral()
            for i in range(self._range[0], self._range[1])
        ]  # Get number of jets by jet pT bins
        self._measN2 = [
            self._f2.Get(
                "{}/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}".format(
                    self._directory2, self._NFIN, i
                )
            ).Integral()
            for i in range(self._range[0], 9)
        ]  # Get number of jets by jet pT bins
        self._measBgN1 = [
            self._f.Get(
                "{}/BgTrkNumberBin/BgTrkNumberBinNFin{:02d}JetPt{:02d}".format(
                    self._directory1, self._NFIN, i
                )
            ).Integral()
            for i in range(self._range[0], self._range[1])
        ]  # Get number of background jets
        self._measBgN2 = [
            self._f2.Get(
                "{}/BgTrkNumberBin/BgTrkNumberBinNFin{:02d}JetPt{:02d}".format(
                    self._directory2, self._NFIN, i
                )
            ).Integral()
            for i in range(self._range[0], 9)
        ]  # Get number of background jets
        self._measRndmBgN1 = [
            self._f.Get(
                "{}/hNumber/hNumberNFin{:02d}".format(self._directory1, self._NFIN)
            ).GetBinContent(7 + i)
            for i in range(self._range[0], self._range[1])
        ]  # Get number of background jets
        self._measRndmBgN2 = [
            self._f2.Get(
                "{}/hNumber/hNumberNFin{:02d}".format(self._directory2, self._NFIN)
            ).GetBinContent(7 + i)
            for i in range(self._range[0], 9)
        ]  # Get number of background jets
        print("Range from {} to {}".format(self._range[0], self._range[1]))
        print(self._measN1)
        print(self._measN2)
        self._measN = [
            self._measN1[i]
            if (i < self._range[1] - self._range[0])
            else self._measN2[i]
            for i in range(0, 9 - self._range[0])
        ]
        self._measBgN = [
            self._measBgN1[i]
            if (i < self._range[1] - self._range[0])
            else self._measBgN2[i]
            for i in range(0, 9 - self._range[0])
        ]
        self._measRndmBgN = [
            self._measRndmBgN1[i]
            if (i < self._range[1] - self._range[0])
            else self._measRndmBgN2[i]
            for i in range(0, 9 - self._range[0])
        ]
        print("MeasN1: {} MeasN2: {}".format(self._measN1, self._measN2))
        print("MeasN: {} ".format(self._measN))

    def getGraphs(self, **kwargs):
        return self.getGraphs2(self._range[0], 8, **kwargs)

    def getJetPt(self, **kwargs):
        try:
            if "dir" in kwargs:
                print(
                    "{0[dir]}/JetPt/JetPtNFin{0[NFin]:02d}".format(
                        {"dir": kwargs["dir"], "NFin": self._NFIN}
                    )
                )
                hist = self._f.Get(
                    "{0[dir]}/JetPt/JetPtNFin{0[NFin]:02d}".format(
                        {"dir": kwargs["dir"], "NFin": self._NFIN}
                    )
                )
            else:
                print(
                    "{0[dir]}/JetPt/JetPtNFin{0[NFin]:02d}".format(
                        {"dir": self._directory2, "NFin": self._NFIN}
                    )
                )
                hist = self._f2.Get(
                    "{0[dir]}/JetPt/JetPtNFin{0[NFin]:02d}".format(
                        {"dir": self._directory2, "NFin": self._NFIN}
                    )
                )
        except rootpy.io.file.DoesNotExist:
            return None
        hist.Scale(1.0, "width")
        return hist

    def getHist(self, name, **kwargs):
        """Return a set of histograms

            Args:
          self: Pointer to the object
          name: Name of histograms
          kwargs: list of keyword arguments
        """
        extra = kwargs.get("extra", "")

        hist1 = [
            self._f.Get(
                "{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}{0[extra]}".format(
                    {
                        "dir": self._directory1,
                        "histname": name,
                        "NFin": self._NFIN,
                        "pT": i,
                        "extra": extra,
                    }
                )
            ).Clone()
            for i in range(self._range[0], self._range[1])
        ]  # Get jT histograms from file an array
        hist2 = [
            self._f2.Get(
                "{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}{0[extra]}".format(
                    {
                        "dir": self._directory2,
                        "histname": name,
                        "NFin": self._NFIN,
                        "pT": i,
                        "extra": extra,
                    }
                )
            ).Clone()
            for i in range(self._range[0], 9)
        ]  # Get jT histograms from file an array
        hist = [
            hist1[i] if (i < self._range[1] - self._range[0]) else hist2[i]
            for i in range(0, 9 - self._range[0])
        ]
        # print('{0[dir]}/{0[histname]}/{0[histname]}NFin{0[NFin]:02d}JetPt{0[pT]:02d}'.format({'dir':self._directory, 'histname':name,'NFin':self._NFIN,'pT':1}))
        jetPt = parse_jet_pt_bins(hist)

        for h, N, bgN, rndmbgN in zip(
            hist, self._measN, self._measBgN, self._measRndmBgN
        ):
            h.Sumw2()
            # print("Rebinning {} by {} in set {} that has {} bins".format(h.GetTitle(), self._rebin, self._name, h.GetNbinsX()))
            h.Rebin(self._rebin)
            if self.properties.get("isWeight", False):
                h.SetLineColor(self.properties.get("color", 1))
                h.SetMarkerColor(self.properties.get("color", 1))
                h.Scale(1.0, "width")
            else:
                if kwargs.get("isBg", False):
                    h.SetLineColor(self.properties.get("color", 1) + 1)
                    h.SetMarkerColor(self.properties.get("color", 1) + 1)
                    h.Scale(1.0 / bgN, "width")
                    print("{} is bg".format(name))
                elif kwargs.get("isRndmBg", False):
                    print("Is random background")
                    h.SetLineColor(self.properties.get("color", 1) + 2)
                    h.SetMarkerColor(self.properties.get("color", 1) + 2)
                    h.Scale(1.0 / rndmbgN, "width")
                    print("Scale by {}".format(rndmbgN))
                else:
                    h.SetLineColor(self.properties.get("color", 1))
                    h.SetMarkerColor(self.properties.get("color", 1))
                    h.Scale(1.0 / N, "width")

            h.SetMarkerStyle(self.properties.get("style", 24))
            h.SetMarkerSize(0.5)
            h.SetLineColor(1)

        if kwargs.get("jetpt", False):
            return hist, jetPt
        else:
            return hist


def compareSets(sets, histname):
    """Draws a comparison between sets of a single histogram

      Args:
      sets: List of sets to be used
      histname: Name of histogram to be drawn
    """
    datasets = [dataset.getHist(histname, jetpt=True) for dataset in sets]
    names = [dataset.name() for dataset in sets]
    fig, axs = defs.makegrid(4, 2, xlog=True, ylog=True, d=d)
    axs = axs.reshape(8)
    axs[1].text(
        0.2,
        0.005,
        d["system"] + "\n" + d["jettype"] + "\n" + d["jetalg"] + "\n" + d["trigger"],
        fontsize=7,
    )

    for i, (set, name) in enumerate(zip(datasets, names)):
        for jT, pT, ax, j in zip(set[0][1:], set[1][1:], axs, range(0, 9)):
            rplt.errorbar(
                jT, xerr=False, emptybins=False, axes=ax, label=name, fmt="o"
            )  # Plot jT histogram,
            if i == 0:
                ax.set_xlim([0.1, 20])  # Set x-axis limits
                ax.set_ylim([5e-4, 2e3])  # Set y-axis limits

    axs[0].legend(loc="lower left")
    plt.show()  # Draw figure on screen


def compareSetsWithRatio(sets, histname):
    """Draw a comparison between sets in 4 jet pT bins with ratios
      to first set in list

      Args:
      sets: List os sets to be used
      histname: Name of histogram to be plotted

    """
    datasets = [dataset.getHist(histname, jetpt=True) for dataset in sets]
    names = [dataset.name() for dataset in sets]
    fig, axs = defs.makegrid(4, 2, xlog=True, ylog=True, d=d, shareY=False)
    axs = axs.reshape(8)
    axs[1].text(
        0.2,
        0.005,
        d["system"] + "\n" + d["jettype"] + "\n" + d["jetalg"] + "\n" + d["trigger"],
        fontsize=7,
    )

    for i, (set, name) in enumerate(zip(datasets, names)):
        for jT, pT, ax, j in zip(set[0][1::2], set[1][1::2], axs[0:4], range(0, 5)):
            rplt.errorbar(
                jT, xerr=False, emptybins=False, axes=ax, label=name, fmt="o"
            )  # Plot jT histogram,
            ax.set_xlim([0.1, 20])  # Set x-axis limits
            ax.set_ylim([5e-4, 2e3])  # Set y-axis limits
            ax.text(
                0.3,
                1e2,
                r"$p_{{T,\mathrm{{jet}}}}$:"
                "\n"
                r" {:02d}-{:02d} GeV".format(pT[0], pT[1]),
            )

    ratiosets = []

    for set in datasets:
        ratios = []
        for s, divider in zip(set[0], datasets[0][0]):
            h = s.Clone()
            h.Divide(divider)
            ratios.append(h)
        ratios1 = (ratios, set[1])
        ratiosets.append(ratios1)

    axs[4].set_ylabel(
        "Ratio", fontsize=9
    )  # Add y-axis labels to left- and righmost subfigures
    axs[-1].set_ylabel("Ratio", fontsize=9)

    for i, (set, name) in enumerate(zip(ratiosets[1:], names[1:])):
        for ratio, pT, ax, j in zip(set[0][1::2], set[1][1::2], axs[4:9], range(0, 5)):
            rplt.errorbar(
                ratio, xerr=False, emptybins=False, axes=ax, label=name, fmt="o"
            )  # Plot ratio histogram,
            # if(i == 0):
            ax.set_xlim([0.1, 20])  # Set x-axis limits
            ax.set_ylim([0.1, 1.8])  # Set y-axis limits
            ax.set_yscale("linear")
    axs[0].legend(loc="lower left")

    return sets


def compareHistsWithRatio(dataset, histnames, labels, step=2, start=1, extras=None):
    """Draw a comparison between sets in 4 jet pT bins with ratios to first set in list

      Args:
      sets: List of sets to be used
      histname: Name of histogram to be plotted

    """
    print(extras)
    if extras is None:
        print("Extras is None")
        hists = [dataset.getHist(histname, jetpt=True) for histname in histnames]
    else:
        print("Extras is ")
        print(extras)
        hists = [
            dataset.getHist(histname, jetpt=True, extra=extra)
            for histname, extra in zip(histnames, extras)
        ]
    fig, axs = defs.makegrid(4, 2, xlog=True, ylog=True, d=d, shareY=False)
    axs = axs.reshape(8)
    axs[1].text(
        0.2,
        0.005,
        d["system"] + "\n" + d["jettype"] + "\n" + d["jetalg"] + "\n" + d["trigger"],
        fontsize=7,
    )
    colors = [1, 2, 3, 4, 5]
    for i, (hist, name) in enumerate(zip(hists, labels)):
        for jT, pT, ax, j in zip(
            hist[0][start::step], hist[1][start::step], axs[0:4], range(0, 5)
        ):
            jT.SetMarkerColor(colors[i])
            rplt.errorbar(
                jT, xerr=False, emptybins=False, axes=ax, label=name, fmt="o"
            )  # Plot jT histogram,
            ax.set_xlim([0.1, 20])  # Set x-axis limits
            ax.set_ylim([5e-4, 2e3])  # Set y-axis limits
            ax.text(
                0.3,
                1e2,
                r"$p_{{T,\mathrm{{jet}}}}$:"
                "\n"
                r" {:02d}-{:02d} GeV".format(pT[0], pT[1]),
            )

    ratiosets = []

    for hist in hists:  # loop over histograms to be compared
        ratios = []
        # divider = hists[0][0]
        for s, divider in zip(
            hist[0], hists[0][0]
        ):  # hist[0] is list of histograms, (hist is a tuple with histograms + jetpT bins)
            h = s.Clone()
            h.Divide(divider)
            ratios.append(h)
        ratios1 = (ratios, hist[1])
        ratiosets.append(ratios1)

    axs[4].set_ylabel(
        "Ratio", fontsize=9
    )  # Add y-axis labels to left- and righmost subfigures
    axs[-1].set_ylabel("Ratio", fontsize=9)

    for i, set in enumerate(ratiosets[1:]):
        for ratio, pT, ax, j in zip(
            set[0][1::step], set[1][1::step], axs[4:9], range(0, 5)
        ):
            rplt.errorbar(
                ratio, xerr=False, emptybins=False, axes=ax, fmt="o"
            )  # Plot ratio histogram,
            # if(i == 0):
            ax.set_xlim([0.1, 20])  # Set x-axis limits
            ax.set_ylim([0.1, 5.8])  # Set y-axis limits
            ax.set_yscale("linear")
    axs[0].legend(loc="lower left")

    return hists


def JtWithBackgroundRatio(dataset, histname, bgname, rndmBgname="", **kwargs):
    """
    Draws a 4x2 grid with jT and Bg on top row for 4 jet pT bins and ratio between these on bottom row

      Args:
      dataset: Dataset to be used
      histname: Name of histogram to be used
      bgname: Name of background histogram to be used
    """
    jtHistos = dataset.getHist(histname, jetpt=True)
    bgHistos = dataset.getHist(bgname, jetpt=True, isBg=True)

    if kwargs.get("includeRandom", False):
        bgRndmHisto = dataset.getHist(rndmBgname, jetpt=True, isRndmBg=True)
    start = kwargs.get("start", 0)
    step = kwargs.get("step", 2)
    fig, axs = defs.makegrid(4, 2, xlog=True, ylog=True, d=d, shareY=False)
    axs = axs.reshape(8)
    axs[1].text(
        1.2,
        10,
        d["system"]
        + "\n"
        + d["jettype"]
        + "\n"
        + d["jetalg"]
        + "\n"
        + d["trigger"]
        + "\n"
        + dataset.name(),
        fontsize=7,
    )

    for jT, bgJt, rndmbg, pT, ax, j in zip(
        jtHistos[0][start::step],
        bgHistos[0][start::step],
        bgRndmHisto[0][start::step],
        jtHistos[1][start::step],
        axs[0:4],
        range(0, 5),
    ):
        rplt.errorbar(
            jT, xerr=False, emptybins=False, axes=ax, label="Inclusive Jt", fmt="o"
        )  # Plot jT histogram,
        rplt.errorbar(
            bgJt, xerr=False, emptybins=False, axes=ax, label="Perp. cone Bg", fmt="o"
        )  # Plot jT histogram,
        if kwargs.get("includeRandom", False):
            rplt.errorbar(
                rndmbg, xerr=False, emptybins=False, axes=ax, label="Random Bg", fmt="o"
            )  # Plot jT histogram,
        ax.set_xlim([0.1, 20])  # Set x-axis limits
        ax.set_ylim([5e-5, 4e3])  # Set y-axis limits
        ax.text(
            0.15,
            5e-4,
            r"$p_{{T,\mathrm{{jet}}}}$:"
            "\n"
            r" {:02d}-{:02d} GeV".format(pT[0], pT[1]),
        )

    ratios = []
    for jT, bg in zip(jtHistos[0], bgHistos[0]):
        h = bg.Clone()
        h.Divide(jT)
        ratios.append(h)
    if kwargs.get("includeRandom", False):
        ratios2 = []
        for jT, bg in zip(jtHistos[0], bgRndmHisto[0]):
            h = bg.Clone()
            h.Divide(jT)
            ratios2.append(h)

    axs[4].set_ylabel(
        "Ratio \n Background/Inclusive", fontsize=9
    )  # Add y-axis labels to left- and righmost subfigures
    axs[-1].set_ylabel("Ratio \n Background/Inclusive", fontsize=9)
    if kwargs.get("includeRandom", False):
        for ratio, pT, ax, j in zip(
            ratios2[1::2], jtHistos[1][1::2], axs[4:9], range(0, 5)
        ):
            rplt.errorbar(
                ratio, xerr=False, emptybins=False, axes=ax, label="Ratio", fmt="o"
            )  # Plot ratio histogram,
    for ratio, pT, ax, j in zip(ratios[1::2], jtHistos[1][1::2], axs[4:9], range(0, 5)):
        rplt.errorbar(
            ratio, xerr=False, emptybins=False, axes=ax, label="Ratio", fmt="o"
        )  # Plot ratio histogram,
        # if(i == 0):
        ax.set_xlim([0.1, 20])  # Set x-axis limits
        ax.set_ylim([0, 2.2])  # Set y-axis limits
        ax.set_yscale("linear")
    axs[0].legend(loc="upper right")


def JtWithBackgroundRatioAll(dataset, histname, bgname):
    """Draws a 4x4 grid with jT and Bg on top and third for 8 jet pT bins and ratio between these on second and bottom row

      Args:
      dataset: Dataset to be used
      histname: Name of histogram to be used
      bgname: Name of background histogram to be used
    """
    jtHistos = dataset.getHist(histname, jetpt=True)
    bgHistos = dataset.getHist(bgname, jetpt=True, isBg=True)
    fig, axs = defs.makegrid(4, 4, xlog=True, ylog=True, d=d, shareY=False)
    axs = axs.reshape(16)
    axs[1].text(
        0.02,
        0.005,
        d["system"] + "\n" + d["jettype"] + "\n" + d["jetalg"] + "\n" + d["trigger"],
        fontsize=7,
    )
    # axs[1].text(0.02, 0.005, d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'] + '\n' + dataset.name(), fontsize = 7)
    print(type(axs))
    for jT, bgJt, pT, ax, j in zip(
        jtHistos[0][1:],
        bgHistos[0][1:],
        jtHistos[1][1:],
        axs[0:4].tolist() + axs[8:12].tolist(),
        range(0, 8),
    ):
        rplt.errorbar(
            jT, xerr=False, emptybins=False, axes=ax, label="Jt", fmt="o"
        )  # Plot jT histogram,
        rplt.errorbar(
            bgJt, xerr=False, emptybins=False, axes=ax, label="BgJt", fmt="o"
        )  # Plot jT histogram,
        ax.set_xlim([0.01, 20])  # Set x-axis limits
        ax.set_ylim([5e-4, 2e3])  # Set y-axis limits
        ax.text(
            0.3,
            1e2,
            r"$p_{{T,\mathrm{{jet}}}}$:"
            "\n"
            r" {:02d}-{:02d} GeV".format(pT[0], pT[1]),
        )

    ratios = []
    for jT, bg in zip(jtHistos[0], bgHistos[0]):
        h = bg.Clone()
        h.Divide(jT)
        ratios.append(h)

    axs[4].set_ylabel(
        "Ratio \n Inclusive/Background", fontsize=9
    )  # Add y-axis labels to left- and righmost subfigures
    axs[-1].set_ylabel("Ratio \n Inclusive/Background", fontsize=9)
    for ratio, pT, ax, j in zip(
        ratios[1:],
        jtHistos[1][1:],
        axs[4:8].tolist() + axs[12:16].tolist(),
        range(0, 8),
    ):
        rplt.errorbar(
            ratio, xerr=False, emptybins=False, axes=ax, label="Ratio", fmt="o"
        )  # Plot ratio histogram,
        # if(i == 0):
        ax.set_xlim([0.01, 20])  # Set x-axis limits
        ax.set_ylim([0, 2.2])  # Set y-axis limits
        ax.set_yscale("linear")
    axs[0].legend(loc="lower left")


def fitJtHisto(histo, method, cut, iJet, iFinder, title="", drawFit=False):
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
    if drawFit:
        defs.drawFit(histo, gaussfit3, iJet, iFinder, title)
        # chi2[iF][ij] = gaussfit3.GetChisquare()
        # chi2[iF][ij] = getChi2(histo, gaussfit3, 0.1, end, 0)
        # chi2dof[iF][ij] = chi2[iF][ij]/(bins[iF][ij]-5)
        # chi2n[iF][ij] = chi2[iF][ij]/bins[iF][ij]
    B2 = gaussfit3.GetParameter(0)
    B1 = gaussfit3.GetParameter(2)
    B2e = gaussfit3.GetParError(0)
    B1e = gaussfit3.GetParError(2)
    gaussigma = math.sqrt(2) * B1
    gausyield = B2 * B1 / math.sqrt(2 * math.pi)
    gausyielde = math.sqrt((B2 * B2 * B1e * B1e + B1 * B1 * B2e * B2e) / (2 * math.pi))
    gaussigmae = math.sqrt(2) * B1e

    constant = gaussfit3.GetParameter(3)  # B3
    constante = gaussfit3.GetParError(3)
    alpha = gaussfit3.GetParameter(5)  # B4
    alphae = gaussfit3.GetParError(5)
    if method == "alt":
        peak = gaussfit3.GetParameter(4)
        peake = gaussfit3.GetParError(4)
        beta = peak * (alpha + 1)  # B5
        betae = math.sqrt(pow((alpha + 1) * peake, 2) + pow(peak * alphae, 2))
    else:
        beta = gaussfit3.GetParameter(4)  # B5
        betae = gaussfit3.GetParError(4)
        peak = beta / (alpha + 1)
        peake = math.sqrt(
            pow(betae / (alpha + 1), 2) + pow(alphae * beta / pow(alpha + 1, 2), 2)
        )

    gammaYield = constant * beta / (alpha - 1)
    gammaRMS = beta / math.sqrt((alpha - 2) * (alpha - 3))
    gammaYielde = math.sqrt(
        pow(beta * constante / (alpha - 1), 2)
        + pow(constant * beta * alphae / pow(alpha - 1, 2), 2)
        + pow(constant * betae / (alpha - 1), 2)
    )
    gammaRMSe = math.sqrt(
        pow(
            (5 - 2 * alpha) * beta * alphae / pow(2 * ((alpha - 2) * (alpha - 3)), 1.5),
            2,
        )
        + pow(betae / math.sqrt((alpha - 2) * (alpha - 3)), 2)
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


def main():
    print("Number of arguments: ", len(sys.argv), "arguments.")
    print("Argument list:", str(sys.argv))
    filename = sys.argv[1]
    print("Input file: ")
    print(filename)
    MB_FullJets_R04 = dataset(
        "MBFullR04",
        NFIN=0,
        filename=filename,
        directory="AliJJetJtTask/AliJJetJtHistManager",
        color=1,
        style=24,
        rebin=5,
    )
    # MB_FullJets_R05 = dataset("MBFullR05", NFIN=1, filename=filename, directory='AliJJetJtTask/AliJJetJtHistManager', color=2, style=24, rebin=5)
    # MB_ChargedJets_R04 = dataset("MBChargedR04", NFIN=2, filename=filename, directory='AliJJetJtTask/AliJJetJtHistManager', color=1, style=24, rebin=5)
    # MB_ChargedJets_R05 = dataset("MBChargedR05", NFIN=3, filename=filename, directory='AliJJetJtTask/AliJJetJtHistManager', color=2, style=24, rebin=5, range = 8)
    # triggered_FullJets_R04 = dataset("TriggeredFullR04", NFIN=0, filename=filename, directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager', color=2, style=24, rebin=5)
    # triggered_FullJets_R05 = dataset("TriggeredFullR05", NFIN=1, filename=filename, directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager', color=1, style=24, rebin=5)
    # triggered_ChargedJets_R04 = dataset("TriggeredChargedR04", NFIN=2, filename=filename, directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager', color=3, style=24, rebin=5)
    # triggered_ChargedJets_R05 = dataset("TriggeredChargedR05", NFIN=3, filename=filename, directory='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager', color=4, style=24, rebin=5)

    # datasets = (MB_FullJets_R04, MB_FullJets_R05, MB_ChargedJets_R04, MB_ChargedJets_R05, triggered_FullJets_R04, triggered_FullJets_R05, triggered_ChargedJets_R04, triggered_ChargedJets_R05)
    # for set in datasets:
    # set.printStats(format ='latex', what = 'jets')
    #   for set in datasets:
    #     set.printStats(format ='latex', what = 'bg')
    #   for set in datasets:
    #     set.printStats(format ='latex', what = 'bgratio')
    #   Mixed_FullJets_R04 = datasetMixed("FullR04", NFIN=0, range=5, filename=filename, directory='AliJJetJtTask/AliJJetJtHistManager', directory2='AliJJetJtTask_kEMCEJE/AliJJetJtHistManager', color=2, style=24, rebin=5)
    #

    # measJt, jetpt = MB_FullJets_R04.getHist('JetConeJtWeightBin', jetpt = True)

    # compareSetsWithRatio((MB_FullJets_R05, MB_FullJets_R04),'JetConeJtWeightBin')
    MB_FullJets_R04.getHistPta("JetConeJtWithPtCutWeightBinBin")
    # compareSetsWithRatio((triggered_FullJets_R05, triggered_FullJets_R04),'JetConeJtWeightBin')
    # Mixed_FullJets_R04.printStats()
    # MB_FullJets_R04.printStats()
    # compareSets((Mixed_FullJets_R04, MB_FullJets_R04),'JetConeJtWeightBin')

    # JtWithBackgroundRatioAll(Mixed_FullJets_R04, 'JetConeJtWeightBin', 'BgJtWeightBin')


# ===============================================================================
#   fig, axs = defs.makegrid(4, 2, xlog=True, ylog=True, d=d)
#   axs = axs.reshape(8)
#   axs[1].text(0.2, 0.005, d['system'] +'\n'+  d['jettype'] +'\n'+ d['jetalg'] + '\n'+ d['trigger'], fontsize = 7)
#
#   for jT, pT, ax, i in zip (measJt[1:], jetpt[1:], axs, range(0, 9)):
#     rplt.errorbar(jT, xerr=False, emptybins=False, axes=ax, label='jT', fmt='+') # Plot jT histogram,
#     ax.set_xlim([0.1, 20]) # Set x-axis limits
#     ax.set_ylim([5e-4, 2e3]) # Set y-axis limits
#     # if i == 0: #For the first subfigure add a legend to bottom left corner
#     #  ax.legend(loc ='lower left')
#   # plt.savefig("jTJetPtBins", format='pdf') #Save figure
#   axs[0].legend(loc = 'lower left')
#   plt.show() # Draw figure on screen
# ===============================================================================


if __name__ == "__main__":
    main()
