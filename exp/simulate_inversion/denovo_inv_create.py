import sys
import pysam
import random
import argparse
import tempfile
import subprocess
import polars as pl
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from typing import Sequence, TypedDict, get_type_hints


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


def run_mm2_dotplot(fa: str, min_aln_len: int) -> pl.DataFrame:
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

        df_paf = pl.read_csv(
            tfh.name,
            has_header=False,
            separator="\t",
            new_columns=get_type_hints(PAF).keys(),
        )
        return df_paf.filter(
            pl.col("aln_len").ge(min_aln_len)
        )


def draw_dotplot(ax: Axes, df_paf: pl.DataFrame):
    first_paf_row: PAF = df_paf.row(0, named=True)
    tlen = first_paf_row["tlen"]
    # Draw self-identity diagonal
    ax.plot(
        [0, tlen],
        [0, tlen],
        color="black"
    )

    for row in df_paf.iter_rows(named=True):
        qst = row["qst"]
        qend = row["qend"]
        tst = row["tst"]
        tend = row["tend"]
        strand = row["strand"]
        if strand == "+":
            x = [qst, qend]
            y = [tst, tend]
            color = "red"
            label = "duplication"
        else:
            x = [tst, tend]
            y = [qend, qst]
            color = "blue"
            label = "inversion"

        ax.plot(x, y, color=color, label=label)


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

    # breakpoint()
    return pysam.FastxRecord(chrom, sequence=final_seq)


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
        help="Minimum aligned block lengths.",
    )
    ap.add_argument(
        "-o", "--output_prefix", type=str, default="./out", help="Output prefix"
    )
    ap.add_argument("-s", "--seed", type=int, default=None, help="Random seed.")
    args = ap.parse_args()

    output_prefix = args.output_prefix

    fasta = pysam.FastaFile(args.fasta)
    assert len(fasta.references) == 1, f"More than one fasta sequence in {args.fasta}"
    chrom = fasta.references[0]

    # Initial dotplot
    print("Running minimap2 to generate initial dotplot", file=sys.stderr)
    df_initial_paf = run_mm2_dotplot(fasta.filename, args.min_aln_len)
    df_initial_paf.write_csv(f"{output_prefix}_before_inv_dotplot.paf", separator="\t")

    fig, axes = plt.subplots(ncols=2, layout="constrained", figsize=(10, 5))
    axes: Sequence[Axes]
    draw_dotplot(axes[0], df_initial_paf)

    fig.savefig(f"{output_prefix}_inv_dotplot.png", bbox_inches="tight")

    # Find regions to swap and induce in fasta
    random.seed(args.seed)
    rand_row = random.randint(0, df_initial_paf.shape[0] - 1)
    row_initial_paf: PAF = df_initial_paf.row(rand_row, named=True)

    inv_seq = generate_inversion(fasta, chrom, row_initial_paf)

    # Then run mm2 again
    with tempfile.NamedTemporaryFile("wt") as tfh:
        tfh.write(str(inv_seq) + "\n")
        tfh.flush()
        df_inv_paf = run_mm2_dotplot(tfh.name, min_aln_len=0)
        df_inv_paf.write_csv(f"{output_prefix}_after_inv_dotplot.paf", separator="\t")

    # Then plot finally
    draw_dotplot(axes[1], df_inv_paf)

    # same axis length since only inversionss
    uniq_labels_handles = {}
    for ax in axes:
        ax.set_xlim(0, row_initial_paf["tlen"])
        ax.set_ylim(0, row_initial_paf["tlen"])

        # Dedup legend
        handles, labels = ax.get_legend_handles_labels()
        labels_handles = dict(zip(labels, handles))
        uniq_labels_handles = uniq_labels_handles | labels_handles

    fig.legend(
        handles=uniq_labels_handles.values(),
        labels=uniq_labels_handles.keys()
    )
    fig.savefig(f"{output_prefix}_inv_dotplot.png", bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
