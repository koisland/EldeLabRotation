import sys
import math
import pysam
import random
import argparse
import tempfile
import subprocess
import polars as pl
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from typing import Sequence, Generator, TypedDict, get_type_hints


nt_table = bytes.maketrans(b"ACTGactg", b"TGACtgac")


class PAF(TypedDict):
    qname: str
    qlen: int
    qst: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tst: int
    tend: int
    matches: int
    aln_len: int
    mapq: int
    tp: str
    cm: str
    s1: str
    dv: str
    rl: str


def run_mm2_dotplot(fa: str) -> pl.DataFrame:
    """
    AvA self-alignment using minimap2 to generate dotplot.

    # Input
    * fa: fasta file of single chromosome

    # Returns
    * `pl.DataFrame` of PAF output
    """
    # https://github.com/lh3/minimap2/issues/106
    # https://lh3.github.io/minimap2/minimap2.html
    # -P - Retain all chains and don’t attempt to set primary chains.
    # -D - If query sequence name/length are identical to the target name/length, ignore diagonal anchors.
    #      This option also reduces DP-based extension along the diagonal.
    # -k - Minimizer k-mer length [15]
    # -w - Minimizer window size [10].
    #      A minimizer is the smallest k-mer in a window of w consecutive k-mers.
    # -m - Discard chains with chaining score <INT [40].
    #      Chaining score equals the approximate number of matching bases minus a concave gap penalty.
    mm2_cmd = ["minimap2", "-PD", "-k19", "-w19", "-m200", "-t8", fa, fa]
    with tempfile.NamedTemporaryFile("wt") as tfh:
        subprocess.run(mm2_cmd, stdout=tfh)
        tfh.flush()

        return pl.read_csv(
            tfh.name,
            has_header=False,
            separator="\t",
            new_columns=get_type_hints(PAF).keys(),
        )


def get_homology_polygons(
    pt_1: tuple[int, int],
    pt_2: tuple[int, int],
    color: str,
    alpha: float,
) -> Generator[Polygon, None, None]:
    rpt_1 = rotate(pt_1)
    rpt_2 = rotate(pt_2)
    coords_l = [rpt_1, rotate((pt_1[0], pt_1[0])), rotate((pt_2[0], pt_2[0])), rpt_2]
    coords_r = [
        rpt_1,
        rotate((pt_1[1], pt_1[1])),
        rotate((pt_2[1], pt_2[1])),
        rpt_2,
    ]
    polygon_l = Polygon(coords_l, facecolor=color, edgecolor="black", alpha=alpha)
    polygon_r = Polygon(coords_r, facecolor=color, edgecolor="black", alpha=alpha)
    yield polygon_l
    yield polygon_r


# https://stackoverflow.com/a/34374437
def rotate(
    point: tuple[int, int],
    origin: tuple[int, int] = (0, 0),
    angle: float = -math.pi / 4,
) -> tuple[float, float]:
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy


def draw_dotplot(ax: Axes, df_paf: pl.DataFrame, min_aln_len: int = 0):
    # filter by alignment length
    df_paf = df_paf.filter(pl.col("aln_len").ge(min_aln_len))

    first_paf_row: PAF = df_paf.row(0, named=True)
    tlen = first_paf_row["tlen"]

    end, _ = rotate((tlen, tlen))

    polygons = []
    row: PAF
    for row in df_paf.iter_rows(named=True):
        tst, tend = row["tst"], row["tend"]
        qst, qend = row["qst"], row["qend"]
        strand = row["strand"]
        if strand == "+":
            pt_1 = (qst, tst)
            pt_2 = (qend, tend)
            color = "red"
            label = "Forward"
        else:
            pt_1 = (tst, qend)
            pt_2 = (tend, qst)
            color = "blue"
            label = "Reverse"

        polygons.extend(get_homology_polygons(pt_1, pt_2, color=color, alpha=0.1))
        pt_1 = rotate(pt_1)
        pt_2 = rotate(pt_2)

        x, y = zip(pt_1, pt_2)
        ax.plot(x, y, color=color, label=label)

    p = PatchCollection(polygons, match_original=True)
    ax.add_collection(p)
    ax.set_ylim(bottom=0)
    ax.set_xlim(0, end)
    ax.set_yticks([], [])

    # Relabel xaxis ticks since rotated (pi / 2)
    ax.xaxis.set_major_formatter(lambda x, pos: f"{x / (math.pi / 2) / 10_000:.1f}")
    ax.set_xlabel("Relative position (Kbp)")


def calculate_midpt(st: int, end: int) -> int:
    midpt = round((end - st) / 2)
    return midpt + st


def generate_inversion(
    fasta: pysam.FastaFile, chrom: str, row_initial_paf: PAF
) -> pysam.FastxRecord:
    seq_len = row_initial_paf["tlen"]
    # Just invert at midpt of homologous sequence
    midpt_1 = calculate_midpt(row_initial_paf["tst"], row_initial_paf["tend"])
    midpt_2 = calculate_midpt(row_initial_paf["qst"], row_initial_paf["qend"])
    if midpt_2 < midpt_1:
        interm = midpt_2
        midpt_2 = midpt_1
        midpt_1 = interm

    print(f"Inverting the region of {midpt_1}-{midpt_2}", file=sys.stderr)
    prev_seq = fasta.fetch(chrom, 0, midpt_1 + 1)
    seq = fasta.fetch(chrom, midpt_1, midpt_2 + 1)
    inv_seq = str(seq.translate(nt_table))
    next_seq = fasta.fetch(chrom, midpt_2, seq_len + 1)
    final_seq = prev_seq + inv_seq + next_seq

    return pysam.FastxRecord(
        f"inv_{chrom}", sequence=final_seq, comment=f"{midpt_1}-{midpt_2}"
    )


def main():
    ap = argparse.ArgumentParser(description="Induce inversion based on dotplot.")
    ap.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=True,
        help="Input fasta file of single region.",
    )
    ap.add_argument(
        "-l",
        "--min_aln_len",
        type=int,
        default=50_000,
        help="Minimum aligned block length to use for inversion.",
    )
    ap.add_argument(
        "-pl",
        "--plot_min_aln_len",
        type=int,
        default=0,
        help="Minimum aligned block length to use for plotting.",
    )
    ap.add_argument("-s", "--seed", type=int, default=None, help="Random seed.")
    ap.add_argument(
        "-o", "--output_prefix", type=str, default="./out", help="Output prefix"
    )
    args = ap.parse_args()

    output_prefix = args.output_prefix
    min_aln_len = args.min_aln_len
    plot_min_aln_len = args.plot_min_aln_len

    fasta = pysam.FastaFile(args.fasta)
    assert len(fasta.references) == 1, f"More than one fasta sequence in {args.fasta}"
    chrom = fasta.references[0]

    # Initial dotplot
    print("Running minimap2 to generate initial dotplot", file=sys.stderr)
    df_initial_paf = run_mm2_dotplot(fasta.filename)
    df_initial_paf.write_csv(f"{output_prefix}_before_inv_dotplot.paf", separator="\t")

    fig, axes = plt.subplots(ncols=2, layout="constrained", figsize=(10, 1.8))
    axes: Sequence[Axes]
    draw_dotplot(axes[0], df_initial_paf, min_aln_len=plot_min_aln_len)

    fig.savefig(f"{output_prefix}_inv_dotplot.png", bbox_inches="tight")

    # Find regions to swap and induce in fasta
    random.seed(args.seed)
    df_subset_initial_paf = df_initial_paf.filter(pl.col("aln_len").ge(min_aln_len))
    if df_subset_initial_paf.is_empty():
        raise RuntimeError(
            f"No valid self-alignments with {min_aln_len=} and fasta file, {fasta.filename}."
        )

    rand_row = random.randint(0, df_subset_initial_paf.shape[0] - 1)
    row_initial_paf: PAF = df_subset_initial_paf.row(rand_row, named=True)

    inv_seq = generate_inversion(fasta, chrom, row_initial_paf)

    # Draw inverted segment on plot
    midpt_1, midpt_2 = [int(pos) for pos in inv_seq.comment.split("-")]

    # Then run mm2 again
    inv_fa = f"{output_prefix}_after_inv.fa"
    with open(inv_fa, "wt") as fh:
        fh.write(str(inv_seq) + "\n")

    df_inv_paf = run_mm2_dotplot(inv_fa)
    df_inv_paf.write_csv(f"{output_prefix}_after_inv_dotplot.paf", separator="\t")

    # Then plot finally
    draw_dotplot(axes[1], df_inv_paf, min_aln_len=plot_min_aln_len)

    # same axis length since only inversionss
    uniq_labels_handles = {}
    for i, ax in enumerate(axes):
        # Dedup legend
        handles, labels = ax.get_legend_handles_labels()
        labels_handles = dict(zip(labels, handles))
        uniq_labels_handles = uniq_labels_handles | labels_handles

        # Draw inverted region
        if i == 0:
            arr_st = (midpt_2, midpt_2)
            arr_end = (midpt_1, midpt_1)
            ax.set_title("Before")
        else:
            arr_st = (midpt_1, midpt_1)
            arr_end = (midpt_2, midpt_2)
            ax.set_title("After")

        ax.annotate(
            "",
            xy=rotate(arr_st),
            xytext=rotate(arr_end),
            arrowprops=dict(facecolor="orange", shrink=0.05),
            label="Inversion",
            zorder=3,
        )

    fig.legend(
        handles=uniq_labels_handles.values(),
        labels=uniq_labels_handles.keys(),
        loc="center left",
        bbox_to_anchor=(1, 0.5),
    )
    fig.savefig(f"{output_prefix}_inv_dotplot.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
