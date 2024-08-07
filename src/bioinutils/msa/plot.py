import matplotlib
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch


sns.set_context("paper", font_scale=1.0, rc={"lines.linewidth": 1.5})
sns.set_style("whitegrid", {"axes.grid": False, "grid.linestyle": ""})
sns.set_style(
    "ticks",
    {
        "axes.grid": False,
        "grid.linestyle": "",
        "font.family": "sans-serif",
        # 'font.sans-serif': 'Myriad Pro',
        "text.color": "0",
        "xtick.color": "0",
        "ytick.color": "0",
    },
)


def check_imgs_path():
    if not os.path.exists("imgs"):
        os.mkdir("imgs")


def apply_nt_to_numeric_to_row(nrow: list):
    cmapper = {"A": 1, "T": 2, "C": 3, "G": 4, "-": 0}
    return [cmapper[c] if c in cmapper else 0 for c in nrow]


def draw_mini_alignment(
    aln_arr: np.array,
    aln_ids: np.array,
    title: str,
    dpi: int = 200,
    letter_size: int = 5,
    height: int = 5,
    width: int = 10,
    show_accession_ids: bool = False,
    annot: bool = False,
    cmap: bool = False,
    show_legend: bool = False,
) -> None:

    cmapper = [
        "white",
        "cornflowerblue",  # A
        "darkred",  # T
        "gold",  # C
        "seagreen",  # G
    ]

    out_name = title.replace(" ", "_")
    out_name = f"imgs/{out_name}_MSA.png"
    # fontsize = 1500 / dpi

    if "-" not in aln_arr:
        cmapper = cmapper[1:]
    if cmap is False:
        cmap = matplotlib.colors.ListedColormap(cmapper)
    caln_arr = np.apply_along_axis(apply_nt_to_numeric_to_row, 1, aln_arr)
    aln_height, aln_width = caln_arr.shape

    fig, ax = plt.subplots(1, 1, dpi=200, figsize=(width, height))
    ax.imshow(caln_arr, cmap=cmap, aspect="auto", interpolation="nearest")
    # ax.set_yticks(list(range(len(aln_ids))))
    # axs.set_xticks(np.arange(1, xlen, 5))

    # thinner lines for high-N alignments
    lineweight_h = 10 / aln_height
    lineweight_v = 10 / aln_width

    # gridlines simple
    ax.hlines(
        np.arange(-0.5, aln_height),
        -0.5,
        aln_width,
        lw=lineweight_h,
        color="white",
        zorder=100,
    )
    ax.vlines(
        np.arange(-0.5, aln_width),
        -0.5,
        aln_height,
        lw=lineweight_v,
        color="white",
        zorder=100,
    )

    # remove borders
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)

    if show_accession_ids is True:
        ax.set_yticks(range(aln_height))
        ax.set_yticklabels(aln_ids)

    plt.title(title)
    plt.ylabel("Sequences")
    plt.xlabel("Alignment Position (bp)")

    if show_legend is True:
        legend_elements = [
            Patch(facecolor=cmapper[-4], label="A"),
            Patch(facecolor=cmapper[-3], label="T"),
            Patch(facecolor=cmapper[-2], label="G"),
            Patch(facecolor=cmapper[-1], label="C"),
        ]

        ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1), loc="center")

    if annot is True:
        for y, seq in enumerate(aln_arr):
            for x, letter in enumerate(seq):
                ax.text(x, y, letter, ha="center", va="center", size=letter_size)

    plt.draw()
    plt.show()
    plt.savefig(out_name, bbox_inches="tight")


sns.set_context("paper", font_scale=1.0, rc={"lines.linewidth": 1.5})
sns.set_style("whitegrid", {"axes.grid": False, "grid.linestyle": ""})
sns.set_style(
    "ticks",
    {
        "axes.grid": False,
        "grid.linestyle": "",
        "font.family": "sans-serif",
        # 'font.sans-serif': 'Myriad Pro',
        "text.color": "0",
        "xtick.color": "0",
        "ytick.color": "0",
    },
)


def plot_pssm(tdf, title):

    title = title.replace(" ", "_")
    check_imgs_path()
    plt.figure(figsize=(14, 3), dpi=200)

    plt.plot(tdf.T)
    plt.title(f"{title} PSSM")
    plt.xlabel("Position (bp)")
    plt.ylabel("Genomes Coverage (%)")
    plt.ylim(0, max(tdf.T.values))
    plt.savefig(f"imgs/{title}_PSSM_percentage.png", bbox_inches="tight")
    plt.draw()


def plot_consensus_threshold(inds, kmer_num, base_file_name):

    check_imgs_path()
    plt.figure(None, figsize=(10, 4))
    plt.plot(inds, kmer_num, color="black")
    plt.title(f"{base_file_name} kmer consensus threshold")
    plt.xlabel("Threshold %")
    plt.ylabel("# Kmers")
    plt.savefig(f"imgs/Consensus_threshold_{base_file_name}.png", bbox_inches="tight")

    plt.draw()
