import logging
import argparse
import polars as pl
import pyideogram as pyid
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.colors import rgb2hex
from matplotlib.patches import Patch


BAND_COLORS = pyid.BANDCOL | {"none": (1.0, 1.0, 1.0)}
LBL_KWARGS = dict(rotation=0, ha="right", va="center")
CHROM_NAMES = [*(str(i) for i in range(1, 23)), "X", "Y"]
BED9_COLS = ["#chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickstart", "thickEnd", "itemRgb"]

logger = logging.getLogger(__name__)


def minimalize_ax(ax: Axes, *, remove_ticks: bool = False) -> None:
    for spine in ["left", "right", "bottom", "top"]:
        ax.spines[spine].set_visible(False)
    if remove_ticks:
        ax.tick_params(
            axis="both",
            left=False,
            top=False,
            right=False,
            bottom=False,
            labelleft=False,
            labeltop=False,
            labelright=False,
            labelbottom=False,
        )


def create_ideogram(
    infile: str,
    annot: str,
    fai: str,
    cytobands: str,
    output_prefix: str,
) -> int:
    df_calls = pl.read_csv(
        infile,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        new_columns=BED9_COLS,
        truncate_ragged_lines=True,
    ).with_columns(length=pl.col("chromEnd") - pl.col("chromStart"))

    df_annot = pl.read_csv(
        annot,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        new_columns=BED9_COLS,
        truncate_ragged_lines=True,
    )

    df_fai = (
        pl.read_csv(
            fai,
            has_header=False,
            separator="\t",
            columns=[0, 1],
            new_columns=["#chrom", "length"]
        )
        .filter(pl.col("#chrom").str.contains("^chr([0-9XY]+)_RagTag"))
        .with_columns(
            chrom_name=pl.col("#chrom").str.extract("^chr([0-9XY]+)")
        )
        .cast({"chrom_name": pl.Enum(CHROM_NAMES)})
        .drop_nulls()
        .sort(by="chrom_name")
        .with_columns(
            row=pl.col("chrom_name").rle_id(),
            col=pl.col("#chrom").str.extract("hap(1|2)").cast(pl.UInt32) - 1
        )
    )
    nrows = df_fai["row"].max() + 1
    df_cytobands = pl.read_csv(
        cytobands,
        separator="\t",
        has_header=False,
        new_columns=["#chrom", "chromStart", "chromEnd"]
    )
    
    color_key = {
        name: rgb2hex([int(e) / 255.0 for e in itemRgb.split(",")])
        if not itemRgb.startswith("#")
        else itemRgb
        for name, itemRgb in df_calls.select("name", "itemRgb").unique().iter_rows()
    }
    max_length = df_fai["length"].max()

    # chrom, calls, spacer
    base_height_ratios = [1.0, 0.25, 1.0]
    width_ratios = [0.5, 0.5]
    num_tracks = len(base_height_ratios)
    height_ratios = base_height_ratios * nrows

    fig, axes = plt.subplots(
        ncols=len(width_ratios),
        nrows=nrows * num_tracks,
        figsize=(20, nrows * 0.75),
        height_ratios=height_ratios,
        width_ratios=width_ratios,
        layout="constrained"
    )

    for row in df_fai.iter_rows(named=True):
        chrom_name = row["#chrom"]
        df_chrom_calls = df_calls.filter(pl.col("#chrom") == chrom_name)
        df_chrom_annot = df_annot.filter(pl.col("#chrom") == chrom_name)
        chrom_length = row["length"]

        ax_row_idx_chrom_track = row["row"] * len(base_height_ratios)
        ax_row_idx_chrom_annot = ax_row_idx_chrom_track + 1
        ax_row_idx_chrom = ax_row_idx_chrom_track + 2
        ax_col_idx_chrom = row["col"]

        ax_chrom: Axes = axes[ax_row_idx_chrom, ax_col_idx_chrom]
        ax_chrom_annot: Axes = axes[ax_row_idx_chrom_annot, ax_col_idx_chrom]
        ax_chrom_track: Axes = axes[ax_row_idx_chrom_track, ax_col_idx_chrom]

        ax_chrom.xaxis.set_tick_params(which="both", length=0, labelleft=False)
        ax_chrom.yaxis.set_tick_params(which="both", length=0)
  
        ax_chrom.set_xlim(0, max_length)
        ax_chrom_annot.set_xlim(0, max_length)
        ax_chrom_track.set_xlim(0, max_length)
        ax_chrom_annot.set_ylim(0, 1)
        ax_chrom_track.set_ylim(0, 1)
        ax_chrom.set_yticks([], [])
        ax_chrom.set_ylabel(chrom_name, **LBL_KWARGS)

        minimalize_ax(ax_chrom, remove_ticks=True)
        minimalize_ax(ax_chrom_annot, remove_ticks=True)
        minimalize_ax(ax_chrom_track, remove_ticks=True)

        # Annot
        for row in df_chrom_annot.iter_rows(named=True):
            color = rgb2hex([int(e) / 255.0 for e in row["itemRgb"].split(",")])
            ax_chrom_annot.axvspan(
                xmin=row["chromStart"], xmax=row["chromEnd"], color=color, label=row["itemRgb"]
            )

        # Write regions
        for row in df_chrom_calls.iter_rows(named=True):
            color = color_key[row["name"]]
            ax_chrom_track.axvspan(
                xmin=row["chromStart"], xmax=row["chromEnd"], color=color, label=row["name"]
            )

        # draw chrom
        ax_chrom.axvspan(
            xmin=0, xmax=chrom_length, color="#d3d3d3"
        )
        for (_, cst, cend) in df_cytobands.filter(pl.col("#chrom") == chrom_name).iter_rows():
            ax_chrom.axvspan(
                xmin=cst, xmax=cend, color="#8b0000"
            )

    # Segdup colors
    legend_patches = [Patch(facecolor=color, label=lbl) for lbl, color in color_key.items()]
    legend_patches.extend([
        Patch(facecolor=color, label=label)
        for label, color in zip(
            [
                r"<90% similarity",
                r"90-98% similarity",
                r"98-99% similarity",
                r">99% similarity",
            ],
            ["#800080", "#808080", "#ffff00ff", "#ffa500"],
        )
    ])
    # Add legend.
    fig.legend(
        handles=legend_patches,
        bbox_to_anchor=(0.5, -0.01),
        loc="upper center",
        ncol=6,
        frameon=False,
        edgecolor="black",
        handlelength=0.7,
        handleheight=0.7,
    )
    # Reduce white space between haps
    logger.info(f"Saving to {output_prefix}.(pdf|png)")
    fig.savefig(f"{output_prefix}.pdf", bbox_inches="tight", dpi=600)
    fig.savefig(f"{output_prefix}.png", bbox_inches="tight", dpi=600)

    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--infile")
    ap.add_argument("-a", "--annot")
    ap.add_argument("-f", "--fai")
    ap.add_argument("-c", "--cytobands")
    ap.add_argument("-o", "--output_prefix")
    args = ap.parse_args()
    return create_ideogram(
        args.infile,
        args.annot,
        args.fai,
        args.cytobands,
        args.output_prefix,
    )


if __name__ == "__main__":
    raise SystemExit(main())
