from __future__ import annotations
import argparse
import logging
import sys
from pathlib import Path

log = logging.getLogger(__name__)


def build_parser():
    p = argparse.ArgumentParser(
        prog="gsp",
        description="Genome-Specific Primer design for polyploid species.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Mode D: BLAST tab + sequence database + genome-wide filter
  python gsp.py -r blast.tab -d sequences.fas -g /data/blast_dbs/jagger \\
                -b bedtools -m muscle -p primer3_core -o out.csv

  # Mode D without genome-wide filter (homeolog check only)
  python gsp.py -r blast.tab -d sequences.fas -b bedtools -m muscle -p primer3_core -o out.csv

  # Mode A: direct multi-FASTA
  python gsp.py -a sequences.fas -g /data/blast_dbs/jagger \\
                -b bedtools -m muscle -p primer3_core -o out.csv

  # Relax 3' end requirement if zero results (try this before lowering -c)
  python gsp.py -r blast.tab -d sequences.fas -g /data/blast_dbs/jagger \\
                -b bedtools -m muscle -p primer3_core -o out.csv --no-require-3end
        """,
    )

    # input modes (mutually exclusive)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("-r", metavar="BLAST_TAB",
                      help="BLAST tabular results file (use with -d)")
    mode.add_argument("-a", metavar="FASTA",
                      help="Multi-FASTA file for direct alignment (no BLAST needed)")

    p.add_argument("-d", metavar="SEQ_DB",
                   help="Sequence database FASTA (required with -r)")
    p.add_argument("-g", metavar="GENOME_BLAST_DB",
                   help="Pre-built BLAST database prefix for genome-wide specificity "
                        "filter (e.g. /data/blast_dbs/jagger). This is the same "
                        "database the user originally BLASTed their query against. "
                        "Omit to skip genome-wide filtering (homeolog check only).")
    p.add_argument("-b", metavar="BEDTOOLS",  required=True,
                   help="Path to bedtools executable")
    p.add_argument("-m", metavar="MUSCLE",    required=True,
                   help="Path to MUSCLE executable (or wrapper script)")
    p.add_argument("-p", metavar="PRIMER3",   required=True,
                   help="Path to primer3_core executable")
    p.add_argument("--blastn", metavar="BLASTN", default="blastn",
                   help="Path to blastn executable (default: blastn in $PATH). "
                        "Only used with -g.")
    p.add_argument("-o", metavar="OUTPUT",    required=True,
                   help="Output CSV file path")
    p.add_argument("-q", metavar="N", type=int, default=3,
                   help="Expected number of subgenome hits per gene (default: 3)")
    p.add_argument("-c", metavar="N", type=int, default=2,
                   help="Minimum unique positions required per primer pair (default: 2)")
    p.add_argument("-t", metavar="SETTINGS",
                   help="Custom Primer3 settings file")
    p.add_argument("-f", metavar="FLANK", type=int, default=1000,
                   help="Flank bp added each side of hit (default: 1000, mode D only)")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Show debug logging")
    p.add_argument("--no-require-3end", action="store_true",
                   help="Relax 3' end uniqueness requirement. By default at least "
                        "one primer per pair must have a unique 3' end base in the "
                        "MSA. Use this only if you get zero results.")
    p.add_argument("--strict-genome", action="store_true",
                   help="In genome-wide BLAST filter, require exactly 1 hit total "
                        "(no exceptions for target chromosomes). Default is polyploid "
                        "mode where hits on the target chromosomes are expected and "
                        "allowed — only hits on OTHER chromosomes disqualify a primer.")
    return p


def validate_args(args):
    """Check argument combinations and file existence."""
    errors = []

    if args.r and not args.d:
        errors.append("-d (sequence database) is required when using -r (BLAST tab)")

    if args.r and not Path(args.r).exists():
        errors.append(f"BLAST tab file not found: {args.r}")

    if args.d and not Path(args.d).exists():
        errors.append(f"Sequence database not found: {args.d}")

    if args.a and not Path(args.a).exists():
        errors.append(f"FASTA file not found: {args.a}")

    if not Path(args.b).exists():
        errors.append(f"bedtools not found: {args.b}")

    if not Path(args.m).exists():
        errors.append(f"MUSCLE not found: {args.m}")

    if not Path(args.p).exists():
        errors.append(f"primer3_core not found: {args.p}")

    if args.t and not Path(args.t).exists():
        errors.append(f"Primer3 settings file not found: {args.t}")

    # Genome BLAST db: check that at least one index file exists
    if args.g:
        db_path = args.g
        nal_file = Path(f"{db_path}.nal")
        if not nal_file.exists():
            errors.append(
                f"Genome BLAST db not found: {nal_file} "
                f"(expected a .nal alias file alongside the folder)"
            )

    if errors:
        for e in errors:
            print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def run_mode_d(args, settings):
    """Mode D: BLAST tab + sequence database."""
    from blast     import parse_blast_tab
    from alignment import build_msa, load_full_sequences_for_group
    from primers   import design_primers
    from output    import write_csv_multi

    require_3end  = not args.no_require_3end
    blast_db      = args.g or None
    blastn_path   = args.blastn

    log.info("Mode D: BLAST tab + sequence database")
    log.info("Parsing BLAST tab: %s", args.r)

    if blast_db:
        log.info("Genome-wide BLAST filter: ENABLED  db=%s", blast_db)
        mode_label = "polyploid (hits on target chroms allowed)"
        if args.strict_genome:
            mode_label = "strict (exactly 1 hit anywhere)"
        log.info("Genome filter mode: %s", mode_label)
    else:
        log.warning(
            "Genome-wide BLAST filter: DISABLED (no -g provided). "
            "Primers will only be checked against the homeolog sequences. "
            "Add -g <blast_db_prefix> to enable full genome specificity."
        )

    if require_3end:
        log.info("3' end uniqueness in MSA: REQUIRED")
    else:
        log.warning("3' end uniqueness: RELAXED (--no-require-3end)")

    groups = parse_blast_tab(args.r)
    log.info("Found %d query group(s)", len(groups))

    results = []

    for group in groups:
        log.info("")
        log.info("═" * 60)
        log.info("Query: %s  (%d chromosomes)", group.query_id, group.num_chroms)
        log.info("═" * 60)

        if group.num_chroms < args.q:
            log.warning(
                "[%s] Only %d chromosomes, need %d — skipping",
                group.query_id, group.num_chroms, args.q
            )
            continue

        selected = group.select_top_n(args.q)
        log.info("Selected chromosomes: %s", selected.chroms)

        # target_chroms = the chromosomes this query hit.
        # In polyploid mode these are the expected (allowed) BLAST hits.
        # In strict mode we pass None so every hit counts against the primer.
        target_chroms = (
            None if args.strict_genome
            else set(selected.chroms)
        )

        try:
            msa, unique = build_msa(selected, args.d, args.m)

            if not unique:
                log.warning("No unique positions found — skipping Primer3")
                continue

            # Load FULL sequences for homeolog-level in-memory check (FIX 1)
            full_seqs = load_full_sequences_for_group(selected, args.d)
            log.info("Full homeolog sequences loaded: %d", len(full_seqs))

            pairs = design_primers(
                msa, unique, args.p,
                settings         = settings,
                min_unique       = args.c,
                full_seqs        = full_seqs,        # FIX 1: full homeolog seqs
                blast_db         = blast_db,          # NEW: genome-wide filter
                blastn_path      = blastn_path,
                target_chroms    = target_chroms,
                require_3end_unique = require_3end,  # FIX 2: 3' end filter
            )

            results.append((group.query_id, pairs))
            log.info("Result: %d primer pairs for %s",
                     len(pairs), group.query_id)

        except Exception as e:
            log.error("Failed on %s: %s", group.query_id, e)
            if args.verbose:
                import traceback
                traceback.print_exc()
            continue

    if not results:
        log.error("No primer pairs produced for any query.")
        log.error("Troubleshooting tips:")
        log.error("  1. Did genome-wide BLAST filter reject everything?")
        log.error("     → Try --strict-genome OFF (default polyploid mode)")
        log.error("     → Check blastn path with: %s -version", blastn_path)
        log.error("  2. Try --no-require-3end to relax MSA 3' end filter")
        log.error("  3. Try -c 1 to lower minimum unique positions threshold")
        log.error("  4. Check your sequences.fas has the right chromosomes")
        sys.exit(1)

    from output import write_csv_multi
    write_csv_multi(results, args.o)
    log.info("")
    log.info("Output: %s", args.o)

    total = sum(len(p) for _, p in results)
    log.info("Total primer pairs written: %d across %d query group(s)",
             total, len(results))


def run_mode_a(args, settings):
    """Mode A: direct multi-FASTA alignment."""
    from alignment import build_msa, read_fasta
    from primers   import design_primers
    from output    import write_csv

    require_3end = not args.no_require_3end
    blast_db     = args.g or None
    blastn_path  = args.blastn

    log.info("Mode A: direct FASTA alignment")
    log.info("Input: %s", args.a)

    if blast_db:
        log.info("Genome-wide BLAST filter: ENABLED  db=%s", blast_db)
    else:
        log.warning(
            "Genome-wide BLAST filter: DISABLED (no -g provided). "
            "Add -g <blast_db_prefix> for full genome specificity."
        )

    # Validate input
    seqs = read_fasta(args.a)
    if not seqs:
        log.error("No sequences found in %s", args.a)
        sys.exit(1)
    if len(seqs) < 2:
        log.error("At least 2 sequences required, only %d found in %s",
                  len(seqs), args.a)
        sys.exit(1)
    log.info("Loaded %d sequences from %s", len(seqs), Path(args.a).name)

    valid = set('ACGTacgtNn-')
    for name, seq in seqs.items():
        invalid = set(seq) - valid
        if invalid:
            log.warning("Sequence '%s' has non-standard characters: %s",
                        name, invalid)

    try:
        msa, unique = build_msa(None, args.a, args.m)
    except Exception as e:
        log.error("Alignment failed: %s", e)
        sys.exit(1)

    log.info("MSA: %d seqs x %d cols, %d unique positions",
             msa.num_seqs, msa.length, len(unique))

    if not unique:
        log.error("No unique positions found — cannot design specific primers.")
        sys.exit(1)

    # In Mode A the input sequences ARE the full sequences —
    # no windowing was applied so pass them directly as full_seqs.
    # target_chroms is None in Mode A (no BLAST hit context).
    pairs = design_primers(
        msa, unique, args.p,
        settings            = settings,
        min_unique          = args.c,
        full_seqs           = seqs,          # FIX 1: full input sequences
        blast_db            = blast_db,       # NEW: genome-wide filter
        blastn_path         = blastn_path,
        target_chroms       = None,           # strict 1-hit mode for Mode A
        require_3end_unique = require_3end,   # FIX 2
    )

    query_id = Path(args.a).stem
    write_csv(pairs, query_id, args.o)
    log.info("Output: %s  (%d pairs)", args.o, len(pairs))


def main():
    parser = build_parser()
    args   = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level  = level,
        format = "%(message)s",
        stream = sys.stderr,
    )

    validate_args(args)

    from primers import load_primer3_settings
    settings = load_primer3_settings(args.t)

    if args.t:
        log.info("Loaded custom Primer3 settings from %s", args.t)

    if args.r:
        run_mode_d(args, settings)
    else:
        run_mode_a(args, settings)


if __name__ == "__main__":
    main()
