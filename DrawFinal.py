from dataset import dataset, d
import defs
import matplotlib.pyplot as plt
import os
import rootpy.plotting.root2matplotlib as rplt
import sys


def main():
    print("Number of arguments: ", len(sys.argv), "arguments.")
    print("Argument list:", str(sys.argv))
    filename = sys.argv[1]
    separate = int(sys.argv[2])
    print("Input file: ")
    print(filename)
    full_jets_r04 = dataset(
        "Full jets R=0.4",
        NFIN=0,
        filename=filename,
        directory="MBDataUnfolder/BayesSubUnfolding",
        color=2,
        style=24,
        rebin=1,
        isWeight=True,
    )

    # The histograms have already been scaled correctly
    # Thus no further scaling is needed
    signal, jetPt = full_jets_r04.getSubtracted(
        "JetConeJtWeightBin", "BgJtWeightBin", jetpt=True, isWeight=True,
    )

    if not os.path.exists('PythonFigures'):
        os.makedirs('PythonFigures')

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
            label=full_jets_r04.name(),
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
                label=full_jets_r04.name(),
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
