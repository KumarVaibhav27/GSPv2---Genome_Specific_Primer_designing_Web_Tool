from __future__ import annotations
import csv
import logging
from pathlib import Path
from primers import PrimerPair

log = logging.getLogger(__name__)

# Column headers matching original GSP output exactly
HEADER = [
    "index", "seq id", "type", "primer start", "primer length",
    "product size",   # note: original has this typo — kept for compatibility
    "Tm", "GC content", "any", "3'", "diff count in primer",
    "3' end diff", "primer seq"
]


def _yes_no(val: bool) -> str:
    return "Yes" if val else "No"


def write_csv(
    pairs:     list[PrimerPair],
    query_id:  str,
    out_path:  str | Path,
    append:    bool = False,
) -> int:
    """
    Write primer pairs to a tab-separated CSV file.
    Each pair produces two rows: LEFT_PRIMER and RIGHT_PRIMER.
    Returns number of pairs written.

    If append=True, appends to existing file (used when multiple
    query groups write to the same output file).
    """
    out_path = Path(out_path)
    mode     = "a" if append else "w"

    written = 0

    with out_path.open(mode, newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")

        # section header comment
        fh.write(f"# Specific primers for {query_id}\n")

        # column header
        writer.writerow(HEADER)

        for rank, pair in enumerate(pairs, 1):
            # LEFT primer row
            writer.writerow([
                rank,
                pair.seq_name,
                "LEFT_PRIMER",
                pair.left_start,
                pair.left_len,
                pair.product_size,
                pair.left_tm,
                pair.left_gc,
                pair.left_any,
                pair.left_end,
                pair.left_diff_count,
                _yes_no(pair.left_3end_diff),
                pair.left_seq,
            ])

            # RIGHT primer row
            writer.writerow([
                rank,
                pair.seq_name,
                "RIGHT_PRIMER",
                pair.right_start,
                pair.right_len,
                pair.product_size,
                pair.right_tm,
                pair.right_gc,
                pair.right_any,
                pair.right_end,
                pair.right_diff_count,
                _yes_no(pair.right_3end_diff),
                pair.right_seq,
            ])

            written += 1

    log.info("Wrote %d primer pairs to %s", written, out_path)
    return written


def write_csv_multi(
    results:  list[tuple[str, list[PrimerPair]]],
    out_path: str | Path,
) -> None:
    """
    Write multiple query groups to a single output file.
    results = [(query_id, pairs), (query_id, pairs), ...]
    First group creates the file, subsequent groups append.
    """
    out_path = Path(out_path)
    for i, (query_id, pairs) in enumerate(results):
        write_csv(pairs, query_id, out_path, append=(i > 0))
        log.info("[%s] %d pairs written", query_id, len(pairs))


# ── self-test ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if len(sys.argv) < 6:
        print("Usage: python output.py <blast.tab> <sequences.fas> "
              "<muscle> <primer3_core> <out.csv>")
        sys.exit(1)

    blast_tab, fasta, muscle_bin, primer3_bin, out_csv = sys.argv[1:6]

    from blast     import parse_blast_tab
    from alignment import build_msa
    from primers   import design_primers

    groups = parse_blast_tab(blast_tab)

    results = []
    for group in groups:
        if group.num_chroms < 2:
            log.warning("Skipping %s — only %d chromosome(s)",
                        group.query_id, group.num_chroms)
            continue
        try:
            msa, unique = build_msa(group, fasta, muscle_bin)
            pairs       = design_primers(msa, unique, primer3_bin)
            results.append((group.query_id, pairs))
            print(f"[{group.query_id}] {len(pairs)} pairs")
        except Exception as e:
            log.error("Failed on %s: %s", group.query_id, e)
            continue

    write_csv_multi(results, out_csv)

    print(f"\nOutput written to: {out_csv}")
    print("First 5 lines:")
    with open(out_csv) as f:
        for i, line in enumerate(f):
            if i >= 5:
                break
            print(" ", line.rstrip())