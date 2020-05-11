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
    print("Number of arguments: ", len(sys.argv), "arguments.")
    print("Argument list:", str(sys.argv))
    filename = sys.argv[1]
    print("Input file: ")
    print(filename)
    MC_Truth = dataset(
        "MCTruth",
        NFIN=0,
        filename=filename,
        directory="Pythia/True",
        color=1,
        style=24,
        rebin=5,
        range=8,
        isWeight=True,
    )
    MC_Obs = dataset(
        "MCObs",
        NFIN=0,
        filename=filename,
        directory="Pythia/Measured",
        color=2,
        style=24,
        rebin=5,
        range=8,
        isWeight=True,
    )
    MC_Unfolded = dataset(
        "Unfolded Bayes",
        NFIN=0,
        filename=filename,
        directory="Pythia/BayesSubUnfolding",
        color=3,
        style=24,
        rebin=5,
        range=8,
        isWeight=True,
    )
    # MC_SVD = dataset("Unfolded SVD",NFIN = 0, filename = filename, directory = 'Pythia/SVDSubUnfolding',color=4,style=24,rebin=5,range=8,isWeight=True)

    compareSetsWithRatio((MC_Truth, MC_Obs, MC_Unfolded), "JetConeJtWeightBin")
    # plt.savefig("PythonFigures/UnfoldingClosure.pdf",format='pdf') #Save figure
    plt.show()  # Draw figure on screen


if __name__ == "__main__":
    main()
