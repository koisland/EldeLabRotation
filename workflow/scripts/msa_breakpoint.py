import sys
import pymafft
import pyfaidx
import argparse
import polars as pl

def main():
    ap = argparse.ArgumentParser(description="Update call bed with scores")
    ap.add_argument("-i", "--infile", help="Input calls")
    ap.add_argument("-f", "--fasta", help="Assembly")
    ap.add_argument("-s", "--bp_slop", default=5, type=int, help="Base pair slop")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose and print msa.")
    args = ap.parse_args()
    bp_slop = args.bp_slop

    fh_fa = pyfaidx.Fasta(args.fasta)
    df_calls = pl.read_csv(
        args.infile,
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom", "st", "end", "name",
            "score", "strand", "tst", "tend", "item_rgb"
        ]
    )
    for row in df_calls.iter_rows(named=True):
        chrom = row["chrom"]
        st = row["st"]
        end = row["end"]
        try:
            brkpt_1 = fh_fa.get_seq(
                chrom,
                st-bp_slop,
                st+bp_slop+1
            )
            brkpt_2 = fh_fa.get_seq(
                chrom,
                end-bp_slop,
                end+bp_slop+1
            )
        except Exception:
            print(f"OOB for ({st-bp_slop}, {st+bp_slop+1}) or ({end-bp_slop}, {end+bp_slop+1})", file=sys.stderr)
            print(*row.values(), sep="\t")
            continue

        msa = pymafft.align(
            [brkpt_1.seq, brkpt_2.seq],
            strategy="auto"
        )
        (_, msa_1), (_, msa_2) = msa.to_tuples()

        max_running_matches = 0
        running_matches = 0
        for nt_1, nt_2 in zip(msa_1, msa_2, strict=True):
            # Count continuous runs of matches
            if nt_1 != "-" and nt_1 == nt_2:
                running_matches += 1
            else:
                running_matches = 0
            max_running_matches = max(max_running_matches, running_matches)

        row["score"] = max_running_matches
        print(*row.values(), sep="\t")

        if args.verbose:
            print(f"Event: {chrom}:{st}-{end}", file=sys.stderr)
            print(f"Breakpoint 1: {chrom}:{st-bp_slop}-{st+bp_slop+1}", file=sys.stderr)
            print(f"Breakpoint 2: {chrom}:{end-bp_slop}-{end+bp_slop+1}", file=sys.stderr)
            print(f"Score: {max_running_matches}", file=sys.stderr)
            print(f" {msa.to_fasta().replace("\n", "\n ")}", file=sys.stderr)
            print("", file=sys.stderr)

if __name__ == "__main__":
    raise SystemExit(main())
