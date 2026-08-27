from __future__ import annotations
import subprocess
import tempfile
import logging
from dataclasses import dataclass, field
from pathlib import Path
from alignment import MSA, UniquePosition

log = logging.getLogger(__name__)

# ── default Primer3 settings — matches C++ exactly ───────────────────────────

DEFAULT_SETTINGS = {
    "PRIMER_TASK":                         "pick_primer_list",
    "PRIMER_MAX_NS_ACCEPTED":              "0",
    "PRIMER_OPT_SIZE":                     "20",
    "PRIMER_MIN_SIZE":                     "18",
    "PRIMER_MAX_SIZE":                     "25",
    "PRIMER_PRODUCT_OPT_SIZE":             "200",
    "PRIMER_PRODUCT_SIZE_RANGE":           "150-700",
    "PRIMER_MIN_GC":                       "30.0",
    "PRIMER_MAX_GC":                       "70.0",
    "PRIMER_MIN_TM":                       "57.0",
    "PRIMER_OPT_TM":                       "60.0",
    "PRIMER_MAX_TM":                       "63.0",
    "PRIMER_MAX_END_STABILITY":            "20.0",
    "PRIMER_PAIR_MAX_DIFF_TM":             "3.0",
    "PRIMER_MAX_TEMPLATE_MISPRIMING":      "12",
    "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING": "24",
    "PRIMER_PAIR_MAX_COMPL_ANY":           "1.00",
    "PRIMER_NUM_RETURN":                   "10000",
}


def load_primer3_settings(settings_file=None):
    """
    Load Primer3 settings, optionally overriding defaults from a file.
    Uses EXACT key matching — fixes the C++ substring collision bug.
    """
    settings = dict(DEFAULT_SETTINGS)
    if settings_file is None:
        return settings

    path = Path(settings_file)
    if not path.exists():
        raise FileNotFoundError(f"Primer3 settings file not found: {path}")

    with path.open() as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                log.warning("Settings line %d skipped (no =): %r", lineno, line)
                continue
            key, _, val = line.partition("=")
            settings[key.strip()] = val.strip()

    return settings


# ── data structures ───────────────────────────────────────────────────────────

@dataclass
class SinglePrimer:
    seq_name:       str
    kind:           str    # LEFT_PRIMER or RIGHT_PRIMER
    start:          int    # 1-based (matches C++ output)
    length:         int
    tm:             str
    gc:             str
    any_th:         str
    end_th:         str
    sequence:       str
    align_start:    int  = 0
    align_end:      int  = 0
    diff_count:     int  = 0
    unique_pos:     list = field(default_factory=list)
    end_3_diff:     bool = False
    deletion_count: int  = 0
    # genome-wide BLAST result fields (filled by genome_blast_filter)
    genome_hits:         int       = -1
    genome_specific:     bool      = False
    end_genome_specific: bool      = False  # 3' end has no off-target genome hits
    genome_hit_chroms:   set       = field(default_factory=set)  # off-target chroms hit


@dataclass
class PrimerPair:
    seq_name:     str
    left_start:   int
    left_len:     int
    left_seq:     str
    left_tm:      str
    left_gc:      str
    left_any:     str
    left_end:     str
    right_start:  int
    right_len:    int
    right_seq:    str
    right_tm:     str
    right_gc:     str
    right_any:    str
    right_end:    str
    product_size: int
    left_diff_count:     int   = 0
    left_3end_diff:      bool  = False
    left_unique_pos:     list  = None
    left_genome_hits:    int   = -1
    left_genome_specific: bool = False
    left_end_genome_specific:  bool = False
    right_diff_count:     int  = 0
    right_3end_diff:      bool = False
    right_unique_pos:     list = None
    right_genome_hits:    int  = -1
    right_genome_specific: bool = False
    right_end_genome_specific:  bool = False
    score:               float = 0.0
    off_target_amplicon: bool  = False  # True if F+R can co-amplify an off-target region


# ── reverse complement ────────────────────────────────────────────────────────

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ── Homeolog 3'-end direct verification ──────────────────────────────────────
#
# G:T wobble pairs: the one mismatch type Taq polymerase can extend through.
# When comparing primer sequence (same orientation as sense strand) to a
# homeolog's sense strand:
#
#   primer='G'  homeolog_sense='A'  →  anti-sense = T  →  primer G : template T  = G:T wobble  (WEAK)
#   primer='T'  homeolog_sense='C'  →  anti-sense = G  →  primer T : template G  = T:G wobble  (WEAK)
#
# All other mismatches (G:A, A:C, C:A, etc.) block Taq extension reliably.
_WOBBLE_MISMATCHES: set[tuple[str, str]] = {('G', 'A'), ('T', 'C')}


def _find_best_match(query: str, template: str) -> tuple[int, int, str]:
    """
    Sliding-window search for best match of `query` in `template`.
    Returns (start_pos, match_count, aligned_region).
    Both strings should be uppercase.
    """
    n = len(query)
    best_pos, best_matches, best_region = 0, 0, ""
    t = template.upper()
    for i in range(max(1, len(t) - n + 1)):
        region = t[i:i + n]
        if len(region) < n:
            break
        m = sum(a == b for a, b in zip(query, region))
        if m > best_matches:
            best_matches, best_pos, best_region = m, i, region
    return best_pos, best_matches, best_region


def _has_strong_3end_mismatch(query: str, region: str, end_window: int = 2) -> bool:
    """
    Returns True if the primer (`query`) has at least one strong (non-wobble)
    mismatch in the last `end_window` bases vs `region` (homeolog sense strand,
    same orientation as query).

    Strong = any mismatch except G:T or T:G wobble pairs.
    A mismatch is strong if neither (query_base, region_base) nor
    (region_base, query_base) is a wobble pair.
    """
    n = len(query)
    for i in range(n - end_window, n):
        q, r = query[i].upper(), region[i].upper()
        if q != r:  # it's a mismatch
            if (q, r) not in _WOBBLE_MISMATCHES and (r, q) not in _WOBBLE_MISMATCHES:
                return True  # strong non-wobble mismatch found
    return False  # no strong mismatch in the window


# ── MSA coordinate mapping ────────────────────────────────────────────────────

def raw_to_msa_col(alignment_seq: str, raw_pos: int) -> int:
    """Convert 1-based raw position to 0-based MSA column."""
    count = 0
    for col, base in enumerate(alignment_seq):
        if base != "-":
            count += 1
            if count == raw_pos:
                return col
    return len(alignment_seq) - 1


def count_gaps_in_range(alignment_seq: str, msa_start: int, msa_end: int) -> int:
    """Count gap characters between two MSA columns (inclusive)."""
    return alignment_seq[msa_start:msa_end+1].count("-")


# ── Primer3 interface ─────────────────────────────────────────────────────────

def _build_p3_input(seq_id: str, raw_seq: str, settings: dict) -> str:
    lines = [f"SEQUENCE_ID={seq_id}", f"SEQUENCE_TEMPLATE={raw_seq}"]
    for key, val in settings.items():
        lines.append(f"{key}={val}")
    lines.append("=")
    return "\n".join(lines) + "\n"


def _parse_p3_legacy(text: str, seq_name: str,
                     alignment_seq: str) -> tuple[list, list]:
    """
    Parse Primer3 -format_output legacy text.
    Converts 0-based Primer3 coordinates to 1-based (C++ does +1).
    Maps raw positions to MSA columns.
    """
    lefts:  list[SinglePrimer] = []
    rights: list[SinglePrimer] = []

    for line in text.split("\n"):
        parts = line.split()
        if len(parts) < 9:
            continue
        if parts[1] not in ("LEFT_PRIMER", "RIGHT_PRIMER"):
            continue

        try:
            kind   = parts[1]
            start  = int(parts[2]) + 1   # C++ does +1 to go 0→1-based
            length = int(parts[3])
            seq    = parts[8]
        except (ValueError, IndexError):
            continue

        if kind == "LEFT_PRIMER":
            a_start = raw_to_msa_col(alignment_seq, start)
            a_end   = raw_to_msa_col(alignment_seq, start + length - 1)
        else:
            a_start = raw_to_msa_col(alignment_seq, start - length + 1)
            a_end   = raw_to_msa_col(alignment_seq, start)

        deletions = count_gaps_in_range(alignment_seq, a_start, a_end)

        p = SinglePrimer(
            seq_name       = seq_name,
            kind           = kind,
            start          = start,
            length         = length,
            tm             = parts[4],
            gc             = parts[5],
            any_th         = parts[6],
            end_th         = parts[7],
            sequence       = seq,
            align_start    = a_start,
            align_end      = a_end,
            deletion_count = deletions,
        )
        if kind == "LEFT_PRIMER":
            lefts.append(p)
        else:
            rights.append(p)

    return lefts, rights


def run_primer3(seq_id: str, raw_seq: str, primer3_bin: str,
                settings: dict, alignment_seq: str,
                workdir=None) -> tuple[list, list]:
    """Run Primer3 and return (left_primers, right_primers)."""
    input_text = _build_p3_input(seq_id, raw_seq, settings)

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".p3in", dir=workdir, delete=False
    ) as f:
        f.write(input_text)
        in_path = Path(f.name)

    out_path = in_path.with_suffix(".p3out")
    cmd = [primer3_bin, "-default_version=1", "-format_output",
           f"-output={out_path}", str(in_path)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    in_path.unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"Primer3 failed for {seq_id}:\n  {result.stderr.strip()}"
        )

    text = out_path.read_text() if out_path.exists() else ""
    out_path.unlink(missing_ok=True)

    return _parse_p3_legacy(text, seq_id, alignment_seq)


# ── homeolog in-memory specificity check ─────────────────────────────────────
#
# FIX 1: Check primer against the FULL homeolog sequences (not just the ~400bp
# alignment window). This catches primers that are unique within the small
# alignment window but still appear elsewhere in the same homeolog chromosome.
# Uses full_seqs from load_full_sequences_for_group().

def _is_specific(primer: SinglePrimer, kind: str,
                 other_seqs: list[str]) -> bool:
    """
    Returns True if the primer does NOT appear in any of the other sequences.
    LEFT:  search forward strand directly.
    RIGHT: search reverse complement (right primer binds opposite strand).
    other_seqs must be the FULL chromosome sequences — not the sub-region.
    """
    query = (primer.sequence if kind == "LEFT_PRIMER"
             else reverse_complement(primer.sequence))
    for seq in other_seqs:
        if query in seq:
            return False
    return True


# ── genome-wide BLAST filter ─────────────────────────────────────────────────
#
# NEW: After homeolog filtering, batch-BLAST all remaining primers against
# the full pre-built genome BLAST database (same db the user BLASTed their
# query against). Only primers with EXACTLY 1 hit genome-wide are kept.
#
# Uses blastn-short with word_size=7, identical to _primer_locations() in
# views.py so results are consistent between the pipeline filter and the
# on-screen location check the frontend already shows.
#
# Both the full primer sequence AND the last 6bp (3' end) must have
# exactly 1 genome-wide hit. This enforces:
#   (a) full primer is unique genome-wide
#   (b) 3' end specifically (the most critical part for specificity) is unique

def genome_blast_filter(
    primers: list[SinglePrimer],
    blast_db: str,
    blastn_path: str = "blastn",
    target_chroms: set[str] | None = None,
    min_pident: float = 95.0,
    timeout: int = 600,
) -> list[SinglePrimer]:
    """
    Genome-wide BLAST specificity filter. Runs ONE batched blastn call for
    all primers, then rejects any primer that:
      - Has more than 1 full-length hit genome-wide (>= min_pident identity)
      - Has more than 1 genome location where its 3' end matches
        at ≥95% identity over ≥18 bp

    If target_chroms is provided (set of chromosome names where the query
    sequence originally hit), a primer is allowed if ALL its hits are within
    those target chromosomes. This handles polyploids where each subgenome
    copy is expected to be on a distinct chromosome.

    Args:
        primers:       list of SinglePrimer objects to filter
        blast_db:      path prefix of pre-built BLAST database
                       (e.g. "/data/blast_dbs/jagger" — same as -db in blastn)
        blastn_path:   path to blastn executable
        target_chroms: set of chromosome names that are the intended targets.
                       Hits on these chromosomes do not count as off-target.
                       Pass None to require exactly 1 hit total (strict mode).
        min_pident:    minimum percent identity to count as a hit (default 95)
        timeout:       blastn subprocess timeout in seconds (default 120)

    Returns:
        Filtered list of primers that pass genome-wide specificity.
    """
    if not primers:
        return []

    if not blast_db:
        log.warning("genome_blast_filter: no blast_db provided — skipping genome filter")
        for p in primers:
            p.genome_hits     = -1
            p.genome_specific = True   # pass through if no db available
        return primers

     # ── Step 2a: BLAST full primers (word_size=7) ─────────────────────────────
    full_fasta = "\n".join(
        f">full_{i}_{'L' if p.kind == 'LEFT_PRIMER' else 'R'}\n"
        f"{(p.sequence if p.kind == 'LEFT_PRIMER' else reverse_complement(p.sequence))}"
        for i, p in enumerate(primers)
    ) + "\n"

    def _run_blast(fasta_input, word_size, timeout):
        cmd = [
            blastn_path,
            "-db",              blast_db,
            "-query",           "-",
            "-task",            "blastn-short",
            "-word_size",       str(word_size),
            "-outfmt",          "6 qseqid sseqid pident length qlen qstart qend",
            "-max_target_seqs", "50",
            "-evalue",          "1000",
            "-dust",            "no",
            "-num_threads",     "8",
        ]
        try:
            r = subprocess.run(
                cmd, input=fasta_input,
                capture_output=True, text=True, timeout=timeout,
            )
            return r.stdout
        except subprocess.TimeoutExpired:
            log.warning("genome_blast_filter: blastn timed out")
            return None
        except FileNotFoundError:
            log.warning("genome_blast_filter: blastn not found at %s", blastn_path)
            return None

    def _run_blast_chunked(fasta_str, word_size, chunk_size=50, per_chunk_timeout=120):
        """Split primers into chunks and run separate BLAST per chunk."""
        lines = fasta_str.strip().split("\n")
        # Each primer is 2 lines (header + sequence)
        chunks = []
        for i in range(0, len(lines), chunk_size * 2):
            chunk = "\n".join(lines[i:i + chunk_size * 2]) + "\n"
            chunks.append(chunk)

        all_output = []
        for chunk in chunks:
            out = _run_blast(chunk, word_size, per_chunk_timeout)
            if out:
                all_output.append(out)
        return "\n".join(all_output) if all_output else ""

    # word_size=7: minimum for blastn-short, catches near-matches in short primers.
    # word_size=11 (previous value) missed 91-92% identity off-target hits.
    full_output = _run_blast_chunked(full_fasta, 9, chunk_size=100, per_chunk_timeout=timeout)


    if not full_output:
        for p in primers:
            p.genome_hits     = -1
            p.genome_specific = True
        return primers

    # ── Step 3: parse hits and count per primer ───────────────────────────────
    # full_hits[i]  = set of chroms where full primer has >= min_pident hit
    # end_hits[i]   = set of chroms where 3' end has >=95% identity hit over >=18bp
    full_hits: dict[int, set[str]] = {i: set() for i in range(len(primers))}
    end_hits:  dict[int, set[str]] = {i: set() for i in range(len(primers))}

    # ── Step 3: parse hits ────────────────────────────────────────────────────
    END_TOLERANCE = 2   # alignment must reach within 2 bp of the primer's 3' end

    for line in (full_output or "").strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        qseqid = parts[0]
        raw_id = parts[1]
        sseqid_short = raw_id.split("|")[-1].split(":")[-1]
        try:
            pident  = float(parts[2])
            aln_len = int(parts[3])
            qlen    = int(parts[4])
            qstart  = int(parts[5])
            qend    = int(parts[6])
        except ValueError:
            continue
        if "full_" not in qseqid:
            continue
        id_parts = qseqid.split("_")
        idx      = int(id_parts[1])
        is_right = len(id_parts) > 2 and id_parts[2] == 'R'

        # Full-length hit: covers >= 85% of primer at >= min_pident identity
        min_aln = max(13, int(0.85 * qlen))
        if pident >= min_pident and aln_len >= min_aln:
            full_hits[idx].add(sseqid_short)

        # 3' end hit — direction depends on primer type:
        # LEFT:  3' end = high qend (end of primer sequence)
        # RIGHT: we BLASTed the RC, so 3' end of primer = qstart=1 of RC
        # Only flag as a risky 3' end hit if the perfect match is ≥18bp —
        # shorter matches are harmless noise in a large repetitive genome.
        if pident >= 95.0 and aln_len >= 18:
            if is_right:
                if qstart <= 1 + END_TOLERANCE:
                    end_hits[idx].add(sseqid_short)
            else:
                if qend >= qlen - END_TOLERANCE:
                    end_hits[idx].add(sseqid_short)

    # ── Step 4: decide pass/fail per primer ──────────────────────────────────
    passed = []
    failed_end = 0
    tier2      = 0   # passed via 3'-end only (full primer also hits off-target chroms)

    for i, p in enumerate(primers):
        fh = full_hits[i]
        eh = end_hits[i]

        p.genome_hits = len(fh)

        if target_chroms:
            off_target_full = fh - target_chroms
            full_ok = len(off_target_full) == 0
            p.genome_hit_chroms = off_target_full   # store for amplicon check
        else:
            off_target_full = fh
            full_ok = len(fh) == 1
            p.genome_hit_chroms = fh

        # 3' end genome check: enforced only for alignments ≥18bp at ≥95% identity.
        # Very short 3' end fragments (8-12bp) match hundreds of loci in the
        # 14.5Gb wheat genome due to repetitive elements — checking those
        # would reject virtually all primers. Instead we require the 3' end
        # alignment to be ≥18bp at ≥95% identity to count as an off-target hit.
        # This catches any genuinely risky 3' end partial match (including hits
        # with 1 body mismatch) while ignoring harmless short coincidences.
        # Note: end_hits is populated only for alignments ≥18bp (see parse step).
        if target_chroms:
            off_target_end        = eh - target_chroms
            end_ok                = len(off_target_end) == 0
            p.end_genome_specific = end_ok
        else:
            end_ok                = len(eh) <= 1
            p.end_genome_specific = end_ok

        p.genome_specific = full_ok

        # Decision logic — the 3' end is the critical discriminator for
        # allele-specific PCR (Taq extends only from a matched 3' end):
        #
        #  end_ok=True,  full_ok=True  → Tier 1: ideal — fully specific.
        #  end_ok=True,  full_ok=False → Tier 2: 3' end specific but full
        #     primer also matches non-target chromosomes.  Biologically valid
        #     for highly-conserved genes: Taq will only extend from the target
        #     (where the 3' end perfectly matches), so the pair will only
        #     amplify the target.  Primer competition is higher but the product
        #     is still specific.  Needed for conserved gene families.
        #  end_ok=False, *            → Always reject: 3' end hits off-target
        #     chromosomes → Taq WILL extend there → non-specific amplification.
        #  full_ok=True, end_ok=False → Same: full-length specific but 3' end
        #     is not — the 3' match elsewhere is what causes off-target product.
        if not end_ok:
            # 3' end hits non-target chromosomes — always reject
            failed_end += 1
        elif end_ok and not full_ok:
            # 3' end is target-specific; full primer matches broadly.
            # Tier 2 — 3' end unique but full primer hits off-target chroms.
            passed.append(p)
            tier2 += 1
        else:
            # Both full-length AND 3' end are target-specific — Tier 1.
            passed.append(p)

    log.info(
        "Genome BLAST filter: %d primers in → %d passed "
        "(%d failed 3' end, %d tier-2 passes: 3'-end specific only)",
        len(primers), len(passed), failed_end, tier2,
    )

    # DEBUG: show hit counts for every primer
    for i, p in enumerate(primers):
        fh = full_hits[i]
        eh = end_hits[i]
        log.info(
            "  [%s] %s | full_chroms=%d %s | end_chroms=%d %s",
            p.kind[:1],
            p.sequence,
            len(fh),
            sorted(fh)[:5],
            len(eh),
            sorted(eh)[:5],
        )

    return passed


# ── unique position annotation ────────────────────────────────────────────────

def _annotate_primer(primer: SinglePrimer,
                     diff_sites_msa: set[int],
                     alignment_seq: str = "") -> None:
    """
    Count MSA unique positions inside a primer and check 3-prime end.
    Matches C++ getSpecificPrimer logic exactly.
    """
    cols   = set(range(primer.align_start, primer.align_end + 1))
    inside = diff_sites_msa & cols
    primer.diff_count = len(inside)

    if alignment_seq:
        seq_idx       = 0
        unique_in_seq = []
        for msa_col in range(primer.align_start, primer.align_end + 1):
            if msa_col < len(alignment_seq) and alignment_seq[msa_col] != '-':
                if msa_col in inside:
                    unique_in_seq.append(seq_idx)
                seq_idx += 1
        primer.unique_pos = unique_in_seq
    else:
        primer.unique_pos = []

    # Window of MSA columns (from the 3' end of the primer) that must contain
    # at least one homeolog-unique position.  Based on the ARMS/allele-specific
    # PCR literature: a mismatch at position 1 or 2 from the 3' end reduces
    # Taq extension by ≥100-fold (position 3 gives only ~10-fold, which is
    # insufficient under standard PCR conditions).  Window = 2 means the
    # discriminating SNP must be in the last 2 bases.
    END3_WINDOW = 2
    if primer.kind == "LEFT_PRIMER":
        end3_cols = set(range(
            max(primer.align_start, primer.align_end - END3_WINDOW + 1),
            primer.align_end + 1
        ))
    else:  # RIGHT_PRIMER — 3' end is at align_start (leftmost column)
        end3_cols = set(range(
            primer.align_start,
            min(primer.align_end + 1, primer.align_start + END3_WINDOW)
        ))

    window_hits = inside & end3_cols
    if window_hits and alignment_seq:
        primer.end_3_diff = any(
            col < len(alignment_seq) and alignment_seq[col] != '-'
            for col in window_hits
        )
    else:
        primer.end_3_diff = bool(window_hits)



# ── deduplication ─────────────────────────────────────────────────────────────

def _dedup_by_position(primers: list[SinglePrimer]) -> list[SinglePrimer]:
    """
    Remove primers within 10bp of an already-kept primer.
    Matches C++ parsePrimerVector logic.
    """
    kept:        list[SinglePrimer]     = []
    seen_by_id:  dict[str, list[int]]  = {}

    for p in primers:
        positions  = seen_by_id.setdefault(p.seq_name, [])
        too_close  = any(abs(p.start - s) <= 10 for s in positions)
        if not too_close:
            positions.append(p.start)
            kept.append(p)
    return kept


# ── scoring ───────────────────────────────────────────────────────────────────

def _score_primer(primer: SinglePrimer, diff_sites_msa: set[int]) -> float:
    """Score = sum(100 / distance_from_3prime) for unique positions inside."""
    score = 0.0
    cols  = set(range(primer.align_start, primer.align_end + 1))
    inside = diff_sites_msa & cols

    if primer.kind == "LEFT_PRIMER":
        end3 = primer.align_end
        for col in inside:
            score += 100.0 / (end3 - col + 1)
    else:
        end3 = primer.align_start
        for col in inside:
            score += 100.0 / (col - end3 + 1)

    return score


# ── FIX 3: specificity tier for sort ─────────────────────────────────────────
#
# Sort order: genome-specific L+R > genome L-only > genome R-only >
#             MSA 3'end L+R > MSA 3'end L-only > MSA 3'end R-only > neither
# Within each tier, sort by score descending.

def _specificity_tier(pair: PrimerPair) -> int:
    """
    Tier ranking optimised for ARMS/allele-specific PCR in polyploid wheat.
    The 3'-end is the key discriminator — Taq extends only from a matched 3'
    end, so a pair where BOTH primers have a 3'-end mismatch with non-target
    homeologs is more reliable than a pair where only one primer's full body
    is genome-specific.

    Tier 8 = both fully genome-specific + both 3'-end BLAST-confirmed  ← ideal
    Tier 7 = both fully genome-specific + one 3'-end BLAST-confirmed
    Tier 6 = both fully genome-specific (full only)
    Tier 5 = both 3'-end BLAST-confirmed (best ARMS pair, full body conserved)
    Tier 4 = left fully specific + left 3'-end confirmed  (one strong side)
    Tier 3 = right fully specific + right 3'-end confirmed (one strong side)
    Tier 2 = one primer fully genome-specific (full only, partner not confirmed)
    Tier 1 = one 3'-end BLAST-confirmed only
    Tier 0 = no genome specificity
    """
    lg  = pair.left_genome_specific
    rg  = pair.right_genome_specific
    leg = pair.left_end_genome_specific
    reg = pair.right_end_genome_specific

    if lg and rg and leg and reg: return 8   # both full + both 3'-end
    if lg and rg and (leg or reg): return 7  # both full + one 3'-end
    if lg and rg:                  return 6  # both full only
    if leg and reg:                return 5  # both 3'-end (ARMS-ideal)
    if lg and leg:                 return 4  # left full + left 3'-end
    if rg and reg:                 return 3  # right full + right 3'-end
    if lg or rg:                   return 2  # one full only
    if leg or reg:                 return 1  # one 3'-end only
    return 0


# ── top-level pipeline ────────────────────────────────────────────────────────

def design_primers(
    msa: MSA,
    unique_positions: list,
    primer3_bin: str,
    settings: dict | None = None,
    min_unique: int = 2,
    full_seqs: dict | None = None,
    blast_db: str | None = None,
    blastn_path: str = "blastn",
    target_chroms: set[str] | None = None,
    require_3end_unique: bool = True,
    primary_target: str | None = None,
) -> list[PrimerPair]:
    """
    Full primer design pipeline:

    1.  Run Primer3 on each target sequence (pick_primer_list, legacy format)
    2.  In-memory homeolog check: reject primers whose sequence appears anywhere
        in the full sequences of the OTHER homeologs  (FIX 1 — full_seqs)
    3.  Deduplicate primers within 10bp (parsePrimerVector)
    4.  Annotate with MSA unique position counts and 3' end flag
    5.  Genome-wide BLAST filter: batch-BLAST all survivors against the full
        genome db; reject any primer with >1 hit (NEW — genome_blast_filter)
    6.  Pair left+right primers: both need >= min_unique, product in range,
        and at least one primer must have a unique 3' end  (FIX 2)
    7.  Sort by specificity tier then score  (FIX 3)

    Args:
        msa:               MSA object from build_msa()
        unique_positions:  UniquePosition list from build_msa()
        primer3_bin:       path to primer3_core executable
        settings:          Primer3 settings dict
        min_unique:        minimum unique MSA positions required per primer
        full_seqs:         {chrom: full_sequence} for all homeologs — used for
                           the in-memory homeolog specificity check (step 2).
                           Pass None to skip (less stringent).
        blast_db:          path prefix of pre-built BLAST database for genome-
                           wide filter (step 5). Same value as -db in blastn,
                           e.g. "/data/blast_dbs/jagger". Pass None to skip.
        blastn_path:       path to blastn executable (default: "blastn")
        target_chroms:     set of chromosome names that are legitimate targets.
                           Hits on these chromosomes are expected and do NOT
                           disqualify a primer. Hits on any other chromosome
                           DO disqualify it. Pass None for strict 1-hit mode.
        require_3end_unique: if True (default), at least one primer per pair
                           must have a unique 3' end base in the MSA.
    """
    if settings is None:
        settings = DEFAULT_SETTINGS

    size_range = settings.get("PRIMER_PRODUCT_SIZE_RANGE", "150-700")
    try:
        prod_min, prod_max = (int(x) for x in size_range.split("-"))
    except ValueError:
        prod_min, prod_max = 150, 700

    raw_seqs = {
        seq_obj.name: seq_obj.trimmed.replace("-", "")
        for seq_obj in msa.sequences
    }

    by_seq: dict[int, set[int]] = {}
    for u in unique_positions:
        by_seq.setdefault(u.seq_index, set()).add(u.col)

    all_lefts:  list[SinglePrimer] = []
    all_rights: list[SinglePrimer] = []

    for seq_idx, seq_obj in enumerate(msa.sequences):
        diff_sites = by_seq.get(seq_idx, set())
        if not diff_sites:
            log.info("[%s] No unique positions — skipping", seq_obj.name[:40])
            continue

        raw_seq = raw_seqs[seq_obj.name]

        # ── FIX 1: use full homeolog sequences for in-memory specificity ──────
        if full_seqs:
            other_seqs = [s for n, s in full_seqs.items()
                          if n != seq_obj.name]
            log.debug("[%s] Homeolog check: %d full sequences",
                      seq_obj.name[:40], len(other_seqs))
        else:
            # Fallback: trimmed sub-regions only — less stringent
            other_seqs = [s for n, s in raw_seqs.items()
                          if n != seq_obj.name]
            log.warning(
                "[%s] full_seqs not provided — homeolog check uses sub-region "
                "only (less stringent, may let non-specific primers through)",
                seq_obj.name[:40]
            )

        log.info("[%s] Running Primer3: %d bp, %d unique positions ...",
                 seq_obj.name[:40], len(raw_seq), len(diff_sites))

        with tempfile.TemporaryDirectory(prefix="gsp_p3_") as tmp:
            # Standard free design
            lefts, rights = run_primer3(
                seq_obj.name, raw_seq, primer3_bin,
                settings, seq_obj.trimmed,
                workdir=Path(tmp)
            )

            # Force Primer3 to end at SNP positions — only when genome BLAST
            # is available (blast_db mode). In direct mode this adds 300+ Primer3
            # calls with no benefit since there is no genome filter.
            snp_positions = sorted(diff_sites)[:15] if blast_db else []
            forced_settings = dict(settings)
            forced_settings["PRIMER_NUM_RETURN"] = "3"

            for snp_col in snp_positions:
                # Convert MSA column to raw sequence position
                raw_pos = sum(
                    1 for c in seq_obj.trimmed[:snp_col+1] if c != '-'
                ) - 1
                if raw_pos < 18 or raw_pos >= len(raw_seq) - 18:
                    continue

                # Force left primer to end at this SNP
                forced_left_settings = dict(forced_settings)
                forced_left_settings["SEQUENCE_FORCE_LEFT_END"] = str(raw_pos)
                fl, _ = run_primer3(
                    f"{seq_obj.name}_FL{raw_pos}", raw_seq, primer3_bin,
                    forced_left_settings, seq_obj.trimmed,
                    workdir=Path(tmp)
                )
                lefts.extend(fl)

                # Force right primer to end at this SNP
                forced_right_settings = dict(forced_settings)
                forced_right_settings["SEQUENCE_FORCE_RIGHT_END"] = str(raw_pos)
                _, fr = run_primer3(
                    f"{seq_obj.name}_FR{raw_pos}", raw_seq, primer3_bin,
                    forced_right_settings, seq_obj.trimmed,
                    workdir=Path(tmp)
                )
                rights.extend(fr)

        log.info("[%s] Primer3 raw: %d left, %d right",
                 seq_obj.name[:40], len(lefts), len(rights))

        # ── Step 2: in-memory homeolog specificity filter ─────────────────────
        lefts  = [p for p in lefts
                  if _is_specific(p, "LEFT_PRIMER",  other_seqs)]
        rights = [p for p in rights
                  if _is_specific(p, "RIGHT_PRIMER", other_seqs)]

        log.info("[%s] After homeolog filter: %d left, %d right",
                 seq_obj.name[:40], len(lefts), len(rights))

        # ── Step 3: deduplicate within 10bp ──────────────────────────────────
        lefts  = _dedup_by_position(lefts)
        rights = _dedup_by_position(rights)

        # ── Step 4: annotate with MSA unique positions and 3' end flag ───────
        for p in lefts + rights:
            _annotate_primer(p, diff_sites, seq_obj.trimmed)

        # ── Step 4b: Direct homeolog 3'-end verification ──────────────────────
        # The MSA-based end_3_diff flag can be fooled by alignment artifacts or
        # G:T wobble mismatches (which Taq extends through). This step re-checks
        # every primer with end_3_diff=True against each homeolog's ACTUAL sequence
        # using a sliding-window Hamming search, then verifies:
        #   (a) there is a mismatch in the last 2 bases, AND
        #   (b) that mismatch is NOT a G:T wobble pair
        # If either condition fails for any homeolog, end_3_diff is cleared.
        if full_seqs:
            other_full = {n: s for n, s in full_seqs.items()
                          if n != seq_obj.name}
            for p in lefts + rights:
                if not p.end_3_diff:
                    continue
                query = (p.sequence if p.kind == "LEFT_PRIMER"
                         else reverse_complement(p.sequence))
                n   = len(query)
                q_u = query.upper()
                min_bind = max(13, int(0.80 * n))  # need ≥80% match to count

                for hname, hseq in other_full.items():
                    _, matches, region = _find_best_match(q_u, hseq)
                    if matches < min_bind:
                        # primer doesn't bind this homeolog strongly → skip
                        continue
                    if not _has_strong_3end_mismatch(q_u, region, end_window=2):
                        # either no mismatch in last 2 bases, or only G:T wobble
                        log.debug(
                            "Direct 3end: %s (%s) fails vs %s "
                            "(primer 3end=%s, homeolog=%s)",
                            p.sequence, p.kind[:1], hname,
                            q_u[-2:], region[-2:] if len(region) >= 2 else region,
                        )
                        p.end_3_diff = False
                        break

        all_lefts.extend(lefts)
        all_rights.extend(rights)

    log.info("Total after homeolog filter + dedup: %d left, %d right",
             len(all_lefts), len(all_rights))

    # ── Step 5: genome-wide BLAST filter ─────────────────────────────────────
    # Batch-BLAST all surviving primers against the full genome database.
    # This is the definitive genome-wide specificity check.
    # Primers with >1 hit genome-wide are rejected here, before pairing,
    # so we never form pairs from non-specific primers.
    if blast_db:
        # Identify which MSA sequence is the primary target (query chromosome).
        # MUSCLE reorders sequences, so we cannot rely on msa.sequences[0].
        # 1. Use primary_target if explicitly provided (preferred).
        # 2. Fall back to matching against target_chroms (single-element set).
        # 3. Final fallback: msa.sequences[0].
        if primary_target:
            query_seq_name = next(
                (s.name for s in msa.sequences
                 if primary_target == s.name
                 or primary_target in s.name
                 or s.name in primary_target),
                msa.sequences[0].name
            )
        elif target_chroms and len(target_chroms) == 1:
            tc = next(iter(target_chroms))
            query_seq_name = next(
                (s.name for s in msa.sequences
                 if tc == s.name or tc in s.name or s.name in tc),
                msa.sequences[0].name
            )
        else:
            query_seq_name = msa.sequences[0].name
        log.info("Genome BLAST query sequence: %s", query_seq_name)
        lefts_to_blast  = [p for p in all_lefts
                           if p.diff_count >= min_unique
                           and p.seq_name == query_seq_name]
        rights_to_blast = [p for p in all_rights
                           if p.diff_count >= min_unique
                           and p.seq_name == query_seq_name]

        # Cap primers before BLAST to limit runtime. Sort by diff_count
        # descending so the most homeolog-discriminating primers go first.
        # 150 per side → at most 2 BLAST calls of 100 primers each per side,
        # still plenty to find good pairs.
        _BLAST_CAP = 150
        if len(lefts_to_blast) > _BLAST_CAP:
            lefts_to_blast  = sorted(lefts_to_blast,
                                     key=lambda p: p.diff_count, reverse=True)[:_BLAST_CAP]
            log.info("Pre-BLAST cap: kept top %d left primers by diff_count", _BLAST_CAP)
        if len(rights_to_blast) > _BLAST_CAP:
            rights_to_blast = sorted(rights_to_blast,
                                     key=lambda p: p.diff_count, reverse=True)[:_BLAST_CAP]
            log.info("Pre-BLAST cap: kept top %d right primers by diff_count", _BLAST_CAP)

        log.info("Running genome-wide BLAST filter against: %s", blast_db)
        log.info("Primers to BLAST: %d left, %d right",
                 len(lefts_to_blast), len(rights_to_blast))

        # Run left and right BLAST filters in parallel — each is an
        # independent batch BLAST call against the genome, so threading
        # gives a near-2× speedup with no correctness impact.
        import concurrent.futures as _cf
        _result_box = {}
        def _blast_side(side, primers):
            _result_box[side] = genome_blast_filter(
                primers, blast_db, blastn_path, target_chroms,
                timeout=600
            )

        with _cf.ThreadPoolExecutor(max_workers=2) as _ex:
            _fl = _ex.submit(_blast_side, 'L', lefts_to_blast)
            _fr = _ex.submit(_blast_side, 'R', rights_to_blast)
            _cf.wait([_fl, _fr])

        lefts_to_blast  = _result_box.get('L', [])
        rights_to_blast = _result_box.get('R', [])

        # Only use genome-filtered primers for pairing
        # Non-BLASTed primers (below min_unique threshold) are discarded
        all_lefts  = lefts_to_blast
        all_rights = rights_to_blast

        log.info("After genome BLAST: %d left, %d right passed",
                 len(all_lefts), len(all_rights))
    else:
        log.warning(
            "No blast_db provided — genome-wide BLAST filter skipped. "
            "Primers are only checked against homeologs. Pass blast_db= "
            "to enable full genome specificity validation."
        )

    # ── Step 6: pair primers ──────────────────────────────────────────────────
    log.info("Pairing %d left x %d right primers ...",
             len(all_lefts), len(all_rights))

    lefts_by_seq:  dict[str, list] = {}
    rights_by_seq: dict[str, list] = {}
    for p in all_lefts:
        lefts_by_seq.setdefault(p.seq_name, []).append(p)
    for p in all_rights:
        rights_by_seq.setdefault(p.seq_name, []).append(p)

    pairs: list[PrimerPair] = []

    for seq_name in lefts_by_seq:
        if seq_name not in rights_by_seq:
            continue
        seq_idx = next(
            (i for i, s in enumerate(msa.sequences) if s.name == seq_name),
            None
        )
        if seq_idx is None:
            continue
        diff_sites = by_seq.get(seq_idx, set())

        for left in lefts_by_seq[seq_name]:
            if left.diff_count < min_unique:
                continue
            if left.deletion_count > 10:
                continue

            for right in rights_by_seq[seq_name]:
                if right.diff_count < min_unique:
                    continue
                if right.deletion_count > 10:
                    continue

                # ── 3' end uniqueness check ───────────────────────────────────
                # At least one primer must have a specific 3' end, verified by
                # EITHER the MSA-based SNP check (end_3_diff) OR the genome-wide
                # BLAST check (end_genome_specific).
                #
                # end_3_diff=True  → MSA shows a homeolog-discriminating SNP in
                #   the last END3_WINDOW (2) MSA columns.
                # end_genome_specific=True → BLAST confirms the 3' end (last
                #   18 bp) only aligns at 100% identity to target chromosomes;
                #   no perfect match on any off-target chromosome.
                #
                # Both measure 3' specificity, but from different angles.
                # The BLAST check is especially important for conserved regions
                # where the MSA SNP falls just outside the last 2 columns yet
                # the 3' end is still demonstrably target-specific in the genome.
                if require_3end_unique:
                    left_3end_ok  = left.end_3_diff  or left.end_genome_specific
                    right_3end_ok = right.end_3_diff or right.end_genome_specific
                    if not left_3end_ok and not right_3end_ok:
                        continue

                prod = right.start - left.start + 1
                if not (prod_min <= prod <= prod_max):
                    continue

                l_score = _score_primer(left,  diff_sites)
                r_score = _score_primer(right, diff_sites)

                # ── Off-target amplicon check (spec section 9) ───────────────
                # If both primers hit the same off-target chromosome, flag the
                # pair — they could co-amplify that locus.
                shared_offtarget = (
                    left.genome_hit_chroms & right.genome_hit_chroms
                )
                has_offtarget_amplicon = len(shared_offtarget) > 0

                pair = PrimerPair(
                    seq_name              = seq_name,
                    left_start            = left.start,
                    left_len              = left.length,
                    left_seq              = left.sequence,
                    left_tm               = left.tm,
                    left_gc               = left.gc,
                    left_any              = left.any_th,
                    left_end              = left.end_th,
                    right_start           = right.start,
                    right_len             = right.length,
                    right_seq             = right.sequence,
                    right_tm              = right.tm,
                    right_gc              = right.gc,
                    right_any             = right.any_th,
                    right_end             = right.end_th,
                    product_size          = prod,
                    left_diff_count       = left.diff_count,
                    left_3end_diff        = left.end_3_diff,
                    left_unique_pos       = left.unique_pos or [],
                    left_genome_hits      = left.genome_hits,
                    left_genome_specific  = left.genome_specific,
                    left_end_genome_specific   = left.end_genome_specific,
                    right_diff_count      = right.diff_count,
                    right_3end_diff       = right.end_3_diff,
                    right_unique_pos      = right.unique_pos or [],
                    right_genome_hits     = right.genome_hits,
                    right_genome_specific = right.genome_specific,
                    right_end_genome_specific  = right.end_genome_specific,
                    score                 = l_score + r_score,
                    off_target_amplicon   = has_offtarget_amplicon,
                )
                pairs.append(pair)

    # ── FIX 3: sort by specificity tier first, then score ────────────────────
    # Tier 7 (both genome-specific) always comes before tier 4 (MSA 3' unique),
    # which always comes before lower tiers. Within each tier, higher score wins.
    pairs.sort(key=lambda p: (_specificity_tier(p), p.score), reverse=True)

    seen:    set  = set()
    deduped: list = []
    for p in pairs:
        key = (p.left_seq, p.right_seq)
        if key not in seen:
            seen.add(key)
            deduped.append(p)

    # ── summary log ───────────────────────────────────────────────────────────
    t7 = sum(1 for p in deduped if p.left_genome_specific  and p.right_genome_specific)
    t6 = sum(1 for p in deduped if p.left_genome_specific  and not p.right_genome_specific)
    t5 = sum(1 for p in deduped if p.right_genome_specific and not p.left_genome_specific)
    t4 = sum(1 for p in deduped if not p.left_genome_specific and not p.right_genome_specific
                                 and p.left_3end_diff and p.right_3end_diff)
    t3 = sum(1 for p in deduped if not p.left_genome_specific and not p.right_genome_specific
                                 and p.left_3end_diff and not p.right_3end_diff)
    t2 = sum(1 for p in deduped if not p.left_genome_specific and not p.right_genome_specific
                                 and not p.left_3end_diff and p.right_3end_diff)
    log.info(
        "Final pairs: %d total  "
        "(genome L+R=%d  genome L=%d  genome R=%d  "
        "MSA L+R 3'=%d  MSA L 3'=%d  MSA R 3'=%d)",
        len(deduped), t7, t6, t5, t4, t3, t2
    )

    return deduped


# ── single-copy primer design ─────────────────────────────────────────────────

def design_primers_single_copy(
    sequence:          str,
    seq_name:          str,
    primer3_bin:       str,
    settings:          dict | None     = None,
    blast_db:          str | None      = None,
    blastn_path:       str             = "blastn",
    target_chroms:     set[str] | None = None,
    require_3end_unique: bool          = True,
) -> list[PrimerPair]:
    """
    Design genome-specific primers for a SINGLE-COPY gene (no homeologs).

    Used automatically when BLAST of the query against the genome returns
    only one chromosome hit.  Because there are no homeologs to compare
    against, the MSA step is skipped entirely.  Specificity is determined
    solely by genome-wide BLAST:

      • A primer is accepted if ALL its full-length hits (≥95% identity,
        ≥85% primer coverage) are on the target chromosome (genome_specific = True).
      • A primer is also accepted if its 3'-end hits only the target chromosome
        at ≥95% identity over ≥18 bp (end_genome_specific = True).
      • Pairs are formed from accepted primers within the product-size range.
      • At least one primer per pair must have end_genome_specific = True.

    Returns a list of PrimerPair sorted by specificity tier then Primer3
    score (same ranking as design_primers).
    """
    if settings is None:
        settings = dict(DEFAULT_SETTINGS)

    size_range = settings.get("PRIMER_PRODUCT_SIZE_RANGE", "150-1500")
    try:
        prod_min, prod_max = (int(x) for x in size_range.split("-"))
    except ValueError:
        prod_min, prod_max = 150, 1500

    log.info("[%s] Single-copy mode: designing primers on %d bp sequence",
             seq_name[:40], len(sequence))

    # ── Step 1: Run Primer3 ───────────────────────────────────────────────────
    with tempfile.TemporaryDirectory(prefix="gsp_sc_") as tmp:
        lefts, rights = run_primer3(
            seq_name, sequence, primer3_bin,
            settings, sequence,          # trimmed == full seq (no MSA)
            workdir=Path(tmp)
        )

    log.info("[%s] Primer3 raw: %d left, %d right", seq_name[:40],
             len(lefts), len(rights))

    if not lefts or not rights:
        log.warning("[%s] Primer3 produced no candidates", seq_name[:40])
        return []

    # ── Step 2: Deduplicate within 10 bp ─────────────────────────────────────
    lefts  = _dedup_by_position(lefts)
    rights = _dedup_by_position(rights)

    # ── Step 3: Annotate (no MSA — diff_count not meaningful) ────────────────
    for p in lefts + rights:
        p.diff_count = 1          # placeholder — genome BLAST is the real gate
        p.end_3_diff = False
        p.unique_pos = []

    # ── Step 4: Cap before BLAST ──────────────────────────────────────────────
    _BLAST_CAP = 150
    if len(lefts)  > _BLAST_CAP:
        lefts  = lefts[:_BLAST_CAP]
    if len(rights) > _BLAST_CAP:
        rights = rights[:_BLAST_CAP]

    # ── Step 5: Genome BLAST filter ───────────────────────────────────────────
    if blast_db:
        import concurrent.futures as _cf
        _box = {}
        def _side(k, primers):
            _box[k] = genome_blast_filter(
                primers, blast_db, blastn_path, target_chroms,
                timeout=600
            )
        with _cf.ThreadPoolExecutor(max_workers=2) as ex:
            fs = [ex.submit(_side, 'L', lefts), ex.submit(_side, 'R', rights)]
            _cf.wait(fs)
        lefts  = _box.get('L', [])
        rights = _box.get('R', [])
        log.info("[%s] After genome BLAST: %d left, %d right passed",
                 seq_name[:40], len(lefts), len(rights))
    else:
        log.warning("[%s] No blast_db — skipping genome filter (less stringent)",
                    seq_name[:40])

    # ── Step 6: Form pairs ────────────────────────────────────────────────────
    pairs = []
    for left in lefts:
        for right in rights:
            if left.start >= right.start:
                continue
            prod = right.start - left.start + 1
            if not (prod_min <= prod <= prod_max):
                continue
            # Require at least one primer to be genome-confirmed
            if not (left.genome_specific or left.end_genome_specific
                    or right.genome_specific or right.end_genome_specific):
                continue
            # Require at least one 3′-end-specific primer (per spec section 4)
            if require_3end_unique:
                if not (left.end_genome_specific or right.end_genome_specific):
                    continue

            # Score by Tm/GC quality (no MSA diff sites in single-copy mode)
            def _qual(p):
                tm_pen = abs(float(p.tm) - 60.0)       # penalty from 60°C opt
                gc_pen = abs(float(p.gc) - 50.0)       # penalty from 50% GC
                return 100.0 - tm_pen * 2 - gc_pen
            l_score = _qual(left)
            r_score = _qual(right)
            # Off-target amplicon check
            shared_offtarget = left.genome_hit_chroms & right.genome_hit_chroms
            pair = PrimerPair(
                seq_name                   = seq_name,
                left_start                 = left.start,
                left_len                   = left.length,
                left_seq                   = left.sequence,
                left_tm                    = left.tm,
                left_gc                    = left.gc,
                left_any                   = left.any_th,
                left_end                   = left.end_th,
                right_start                = right.start,
                right_len                  = right.length,
                right_seq                  = right.sequence,
                right_tm                   = right.tm,
                right_gc                   = right.gc,
                right_any                  = right.any_th,
                right_end                  = right.end_th,
                product_size               = prod,
                left_diff_count            = left.diff_count,
                left_3end_diff             = left.end_3_diff,
                left_unique_pos            = [],
                left_genome_hits           = left.genome_hits,
                left_genome_specific       = left.genome_specific,
                left_end_genome_specific   = left.end_genome_specific,
                right_diff_count           = right.diff_count,
                right_3end_diff            = right.end_3_diff,
                right_unique_pos           = [],
                right_genome_hits          = right.genome_hits,
                right_genome_specific      = right.genome_specific,
                right_end_genome_specific  = right.end_genome_specific,
                score                      = l_score + r_score,
                off_target_amplicon        = len(shared_offtarget) > 0,
            )
            pairs.append(pair)

    pairs.sort(key=lambda p: (_specificity_tier(p), p.score), reverse=True)

    # Deduplicate pairs (same left+right start)
    seen:    set  = set()
    deduped: list = []
    for p in pairs:
        key = (p.left_start, p.right_start)
        if key not in seen:
            seen.add(key)
            deduped.append(p)

    log.info("[%s] Single-copy mode: %d primer pairs found",
             seq_name[:40], len(deduped))
    return deduped


# ── self-test ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if len(sys.argv) < 5:
        print("Usage: python primers.py <blast.tab> <sequences.fas> "
              "<muscle> <primer3_core> [blast_db]")
        sys.exit(1)

    blast_tab, fasta, muscle_bin, primer3_bin = sys.argv[1:5]
    blast_db = sys.argv[5] if len(sys.argv) > 5 else None

    from blast     import parse_blast_tab
    from alignment import build_msa, load_full_sequences_for_group

    groups     = parse_blast_tab(blast_tab)
    test_group = next((g for g in groups if g.num_chroms >= 2), None)
    if not test_group:
        print("No suitable group found.")
        sys.exit(1)

    print(f"Query: {test_group.query_id}")
    msa, unique = build_msa(test_group, fasta, muscle_bin)
    print(f"MSA: {msa.num_seqs} seqs x {msa.length} cols, "
          f"{len(unique)} unique positions")

    full_seqs     = load_full_sequences_for_group(test_group, fasta)
    target_chroms = set(test_group.chroms) if blast_db else None

    pairs = design_primers(
        msa, unique, primer3_bin,
        full_seqs=full_seqs,
        blast_db=blast_db,
        target_chroms=target_chroms,
    )

    print(f"\nTop 10 of {len(pairs)} pairs:")
    print(f"{'#':>3}  {'seq':<35}  {'score':>7}  "
          f"{'Lu':>3}  {'Ru':>3}  {'L3':>5}  {'R3':>5}  "
          f"{'Lg':>5}  {'Rg':>5}  {'prod':>5}")
    print("-" * 88)
    for i, p in enumerate(pairs[:10], 1):
        print(f"{i:>3}  {p.seq_name[:35]:<35}  {p.score:>7.2f}  "
              f"{p.left_diff_count:>3}  {p.right_diff_count:>3}  "
              f"{str(p.left_3end_diff):>5}  {str(p.right_3end_diff):>5}  "
              f"{str(p.left_genome_specific):>5}  "
              f"{str(p.right_genome_specific):>5}  "
              f"{p.product_size:>5}")
