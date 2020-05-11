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
    print "Number of arguments: ", len(sys.argv), "arguments."
    print "Argument list:", str(sys.argv)
    filename = sys.argv[1]
    separate = int(sys.argv[2])
    print "Input file: "
    print filename
    FullJets_R04 = dataset(
        "Full jets R=0.4",
        NFIN=0,
        filename=filename,
        directory="data/BayesSubUnfolding",
        color=2,
        style=24,
        rebin=5,
    )
    signal, jetPt = Mixed_FullJets_R04.getSubtracted(
        "JetConeJtWeightBin", "BgJtWeightBin", jetpt=True
    )

    if separate > 0:
        fig = plt.figure(1)
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.text(
            0.2,
            0.0005,
            d["system"] + "\n" + d["jettype"] + "\n" + d["jetalg"] + "\n Jet Cone",
            fontsize=7,
        )
        rplt.errorbar(
            signal[separate],
            xerr=False,
            emptybins=False,
            axes=ax,
            label=Mixed_FullJets_R04.name(),
            fmt="o",
        )  # Plot jT histogram,
        ax.text(
            0.3,
            1e2,
            r"$p_{{T,\mathrm{{jet}}}}$:"
            "\n"
            r" {:02d}-{:02d} GeV".format(jetPt[separate][0], jetPt[separate][1]),
        )
        ax.set_xlim([0.1, 12])
        ax.set_ylim([5e-6, 1.5e3])
        ax.legend(loc="lower left")

        plt.savefig(
            "PythonFigures/UnfoldedFullJetsR04JetConeJtFinalJetPt{0}.pdf".format(
                separate
            ),
            format="pdf",
        )  # Save figure
        plt.show()  # Draw figure on screen

    else:
        fig, axs = defs.makegrid(4, 2, xlog=True, ylog=True, d=d, shareY=True)
        axs = axs.reshape(8)
        axs[1].text(
            0.12,
            0.002,
            d["system"] + "\n" + d["jettype"] + "\n" + d["jetalg"] + "\n Jet Cone",
            fontsize=7,
        )
        for jT, pT, ax, i in zip(signal[1:], jetPt[1:], axs, range(0, 9)):
            rplt.errorbar(
                jT,
                xerr=False,
                emptybins=False,
                axes=ax,
                label=Mixed_FullJets_R04.name(),
                fmt="o",
            )  # Plot jT histogram,

            ax.text(
                0.3,
                1e2,
                r"$p_{{T,\mathrm{{jet}}}}$:"
                "\n"
                r" {:02d}-{:02d} GeV".format(pT[0], pT[1]),
            )

            ax.set_xlim([0.1, 13])  # Set x-axis limits
            ax.set_ylim([5e-4, 2e3])  # Set y-axis limits

        # axs[0].legend(loc = 'lower left')

        plt.savefig(
            "PythonFigures/UnfoldedFullJetsR04JetConeJtFinal.pdf", format="pdf"
        )  # Save figure
        plt.show()  # Draw figure on screen


if __name__ == "__main__":
    main()
