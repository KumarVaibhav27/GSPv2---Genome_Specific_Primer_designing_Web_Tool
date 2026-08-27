from __future__ import annotations
import subprocess
import tempfile
import logging
from dataclasses import dataclass
from pathlib import Path
from blast import HitGroup

log = logging.getLogger(__name__)

# ── data structures ──────────────────────────────────────────────────────────

@dataclass
class AlignedSequence:
    name:    str
    raw:     str
    aligned: str
    trimmed: str

@dataclass
class MSA:
    query_id:  str
    sequences: list

    @property
    def length(self):
        return len(self.sequences[0].trimmed) if self.sequences else 0

    @property
    def num_seqs(self):
        return len(self.sequences)

    def column(self, i):
        return [s.trimmed[i] for s in self.sequences]

@dataclass
class UniquePosition:
    col:         int
    seq_index:   int
    seq_name:    str
    unique_base: str

# ── FASTA helpers ─────────────────────────────────────────────────────────────

def read_fasta(path):
    """
    Read any FASTA file into {header: sequence}.
    Handles both gapped (alignment) and ungapped sequences.
    """
    result, current, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current is not None:
                    result[current] = "".join(parts).upper()
                current = line[1:].strip()
                parts = []
            elif line:
                parts.append(line)
    if current is not None:
        result[current] = "".join(parts).upper()
    return result

# ── BED coordinate calculation (matches C++ generateBedLocationFile) ─────────

def _calc_bed_coords(hits, flank, seq_len):
    """
    Calculate extraction coordinates for a chromosome given its BLAST hits.
    Matches C++ generateBedLocationFile logic exactly:
    - Collect min/max of all hit start/end coordinates
    - Add flank on each side, clamp to [1, seq_len]
    - If span > 30000bp, fall back to first hit +/- flank only
    Returns (start, end) as 1-based inclusive coordinates.
    """
    starts = [h.sstart for h in hits]
    ends   = [h.send   for h in hits]
    all_coords = starts + ends

    min_coord = min(all_coords)
    max_coord = max(all_coords)

    tmp_min = min_coord - flank
    tmp_max = max_coord + flank

    if tmp_max - tmp_min > 30000:
        first_start = hits[0].sstart
        first_end   = hits[0].send
        if first_start < first_end:
            tmp_min = first_start - flank
            tmp_max = first_end   + flank
        else:
            tmp_min = first_end   - flank
            tmp_max = first_start + flank

    tmp_min = max(1, tmp_min)
    tmp_max = min(seq_len, tmp_max)
    return tmp_min, tmp_max


# ── sequence loading ──────────────────────────────────────────────────────────

def load_sequences_for_group(group, fasta_path, min_chroms=2):
    """
    Load sequences from sequences.fas for all chromosomes in the hit group.
    Mode -a: use all sequences in the file.
    Mode -d: extract sub-region around BLAST hit (hit +/- flank).
    Returns {chrom_name: sequence}.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"Sequence FASTA not found: {fasta_path}")

    all_seqs = read_fasta(fasta_path)
    log.debug("FASTA contains %d sequences", len(all_seqs))

    if group is None:
        log.info("Mode A: loaded %d sequences from %s",
                 len(all_seqs), fasta_path.name)
        return all_seqs

    FLANK = 200
    selected = {}
    missing  = []

    for chrom, hits in group.hits_by_chrom.items():
        if chrom not in all_seqs:
            missing.append(chrom)
            log.warning("Chromosome not in FASTA, skipping: %s", chrom)
            continue

        full_seq = all_seqs[chrom]
        seq_len  = len(full_seq)
        bed_start, bed_end = _calc_bed_coords(hits, FLANK, seq_len)
        sub_seq = full_seq[bed_start - 1 : bed_end]

        strand = hits[0].strand
        if strand == "-":
            comp   = str.maketrans("ACGTacgt", "TGCAtgca")
            sub_seq = sub_seq.translate(comp)[::-1]

        log.debug("%s: bed=[%d,%d] strand=%s len=%d",
                  chrom[:40], bed_start, bed_end, strand, len(sub_seq))
        selected[chrom] = sub_seq

    if missing:
        log.warning("[%s] %d/%d chromosomes not found in FASTA: %s",
                    group.query_id, len(missing), group.num_chroms,
                    ", ".join(m[:40] for m in missing))

    if len(selected) < min_chroms:
        raise ValueError(
            f"[{group.query_id}] Only {len(selected)} chromosome(s) found "
            f"in FASTA (need >= {min_chroms}). Missing: {missing}"
        )

    log.info("[%s] Loaded %d/%d chromosome sub-regions",
             group.query_id, len(selected), group.num_chroms)
    return selected


def load_full_sequences_for_group(group, fasta_path):
    """
    Load the COMPLETE unclipped sequences for each chromosome in the group.

    Used for homeolog-level specificity checking (in-memory string search).
    Ensures a primer is checked against the full homeolog sequence, not just
    the ~400bp alignment window — preventing primers that match outside the
    window from slipping through.

    For genome-wide specificity, use genome_blast_filter() instead which
    checks against the full pre-built BLAST database.

    Returns {chrom_name: full_sequence}.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"Sequence FASTA not found: {fasta_path}")

    all_seqs = read_fasta(fasta_path)

    if group is None:
        return all_seqs

    full = {}
    for chrom in group.hits_by_chrom:
        if chrom in all_seqs:
            full[chrom] = all_seqs[chrom]
        else:
            log.warning("Full seq not found for specificity check: %s", chrom)

    return full


# ── MUSCLE alignment ──────────────────────────────────────────────────────────

_MUSCLE_MAX_BP = 10_000   # MUSCLE 5 osxarm64 crashes silently above ~12 kb

def align_sequences(sequences, muscle, workdir=None):
    """
    Write sequences to temp FASTA, run MUSCLE, return {name: aligned_seq}.

    Sequences longer than _MUSCLE_MAX_BP are trimmed to their central window
    before alignment. MUSCLE 5 on Apple Silicon silently exits non-zero for
    inputs above ~12 kb; trimming keeps alignment reliable while still
    covering any primer-product window (amp_max ≤ 1500 bp << _MUSCLE_MAX_BP).
    """
    # Trim oversized sequences to the central _MUSCLE_MAX_BP window
    longest = max(len(s.replace('-', '')) for s in sequences.values())
    if longest > _MUSCLE_MAX_BP:
        log.info(
            "align_sequences: sequences up to %d bp — trimming to central "
            "%d bp to avoid MUSCLE crash on large inputs", longest, _MUSCLE_MAX_BP
        )
        trimmed = {}
        for name, seq in sequences.items():
            raw = seq.replace('-', '')
            if len(raw) > _MUSCLE_MAX_BP:
                mid  = len(raw) // 2
                half = _MUSCLE_MAX_BP // 2
                raw  = raw[mid - half : mid + half]
            trimmed[name] = raw
        sequences = trimmed

    with tempfile.NamedTemporaryFile(
        mode="w", suffix="_input.fas", dir=workdir, delete=False
    ) as inf:
        for name, seq in sequences.items():
            inf.write(f">{name}\n{seq.replace('-','')}\n")
        in_path = Path(inf.name)

    out_path = in_path.with_suffix(".aligned.fas")
    cmd      = [muscle, "-align", str(in_path), "-output", str(out_path)]

    log.debug("MUSCLE: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

    if result.returncode != 0:
        raise RuntimeError(
            f"MUSCLE failed (exit {result.returncode}):\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  stderr: {result.stderr.strip()}"
        )

    aligned = read_fasta(out_path)
    in_path.unlink(missing_ok=True)
    out_path.unlink(missing_ok=True)

    if not aligned:
        raise RuntimeError(
            f"MUSCLE produced no output. "
            f"Input had {len(sequences)} sequences."
        )

    log.debug("MUSCLE aligned %d sequences", len(aligned))
    return aligned

# ── MSA trimming ──────────────────────────────────────────────────────────────

def trim_alignment(aligned):
    """
    Remove leading/trailing columns where fewer than 2 sequences have a
    real base (non-gap character).
    Matches C++ parseAlignment logic exactly.
    Returns {name: trimmed_seq}.
    """
    if not aligned:
        return {}

    names  = list(aligned.keys())
    seqs   = list(aligned.values())
    n      = len(seqs)
    length = len(seqs[0])

    for name, seq in aligned.items():
        if len(seq) != length:
            raise ValueError(
                f"Ragged alignment: {name!r} has {len(seq)} cols, "
                f"expected {length}"
            )

    def at_least_two_real(col):
        gaps = sum(1 for s in seqs if s[col] == "-")
        return gaps <= n - 2

    left = 0
    while left < length and not at_least_two_real(left):
        left += 1

    right = length - 1
    while right >= left and not at_least_two_real(right):
        right -= 1

    trimmed = {name: seqs[i][left:right+1] for i, name in enumerate(names)}
    log.debug("Trimmed: %d -> %d columns", length, right - left + 1)
    return trimmed

# ── unique position detection ─────────────────────────────────────────────────

def find_unique_positions(msa):
    """
    Find columns where exactly one sequence has a unique character.
    Matches C++ identifyDifferentSite exactly.
    Returns list of UniquePosition sorted by column.
    """
    unique_positions = []

    for col in range(msa.length):
        bases = msa.column(col)

        for seq_index, base in enumerate(bases):
            char_count = bases.count(base)
            if char_count == 1:
                unique_positions.append(UniquePosition(
                    col        = col,
                    seq_index  = seq_index,
                    seq_name   = msa.sequences[seq_index].name,
                    unique_base = base,
                ))

    return unique_positions

# ── top-level pipeline ────────────────────────────────────────────────────────

def build_msa(group, fasta_path, muscle, min_chroms=2):
    """
    Full pipeline for one query group:
      load sequences -> align -> trim -> find unique positions
    Returns (MSA, list[UniquePosition]).
    """
    query_id = group.query_id if group else Path(fasta_path).stem

    raw_seqs = load_sequences_for_group(
        group, fasta_path, min_chroms=min_chroms
    )

    log.info("[%s] Aligning %d sequences with MUSCLE ...",
             query_id, len(raw_seqs))
    log.info("[%s] Sequence names: %s", query_id, list(raw_seqs.keys()))

    if len(raw_seqs) < 2:
        raise RuntimeError(
            f"Need at least 2 sequences for alignment, "
            f"got {len(raw_seqs)}: {list(raw_seqs.keys())}"
        )

    lengths        = [len(seq) for seq in raw_seqs.values()]
    has_gaps       = any('-' in seq for seq in raw_seqs.values())
    all_same_length = len(set(lengths)) == 1

    if has_gaps and all_same_length and len(raw_seqs) >= 2:
        log.info("[%s] Pre-aligned sequences detected — skipping MUSCLE", query_id)
        aligned_seqs = dict(raw_seqs)
    else:
        with tempfile.TemporaryDirectory(prefix="gsp_") as tmpdir:
            aligned_seqs = align_sequences(raw_seqs, muscle, workdir=Path(tmpdir))

    trimmed_seqs = trim_alignment(aligned_seqs)

    seq_objects = [
        AlignedSequence(
            name    = name,
            raw     = raw_seqs.get(name, ""),
            aligned = aligned_seqs.get(name, ""),
            trimmed = trimmed_seqs[name],
        )
        for name in trimmed_seqs
    ]

    msa    = MSA(query_id=query_id, sequences=seq_objects)
    unique = find_unique_positions(msa)

    log.info("[%s] MSA: %d seqs x %d cols, %d unique positions",
             query_id, msa.num_seqs, msa.length, len(unique))

    return msa, unique

# ── self-test ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if len(sys.argv) < 4:
        print("Usage: python alignment.py <blast.tab> <sequences.fas> <muscle>")
        print("   or: python alignment.py --direct <sequences.fas> <muscle>")
        sys.exit(1)

    if sys.argv[1] == "--direct":
        fasta_path = sys.argv[2]
        muscle_bin = sys.argv[3]
        print(f"\nMode A (direct FASTA): {fasta_path}")
        msa, unique = build_msa(None, fasta_path, muscle_bin)
    else:
        blast_tab, fasta_path, muscle_bin = sys.argv[1:4]
        from blast import parse_blast_tab
        groups     = parse_blast_tab(blast_tab)
        test_group = next(
            (g for g in groups if g.num_chroms >= 2), None
        )
        if not test_group:
            print("No group with >= 2 chromosomes found.")
            sys.exit(1)
        print(f"\nMode D (blast+db): query={test_group.query_id}")
        msa, unique = build_msa(test_group, fasta_path, muscle_bin)

    print(f"\nMSA:")
    print(f"  sequences : {msa.num_seqs}")
    print(f"  columns   : {msa.length}")
    for s in msa.sequences:
        gaps = s.trimmed.count("-")
        pct  = 100*gaps//len(s.trimmed) if s.trimmed else 0
        print(f"  {s.name[:55]}: {len(s.trimmed)} cols, {gaps} gaps ({pct}%)")

    print(f"\nUnique positions: {len(unique)}")
    for u in unique[:15]:
        print(f"  col {u.col:5d}  seq={u.seq_index}  "
              f"base={u.unique_base}  [{u.seq_name[:45]}]")
    if len(unique) > 15:
        print(f"  ... and {len(unique)-15} more")
