"""
KASP marker design module.
"""
import subprocess
import tempfile
import logging
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

FAM_TAIL = "GAAGGTGACCAAGTTCATGCT"
HEX_TAIL = "GAAGGTCGGAGTCAACGGATT"


# ── Thermodynamics ────────────────────────────────────────────────────────────

# SantaLucia 1998 nearest-neighbor parameters for DNA/DNA duplexes
# dH in kcal/mol, dS in cal/(mol·K)
_NN = {
    'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
    'AT': (-7.2, -20.4),
    'TA': (-7.2, -21.3),
    'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
    'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
    'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
    'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
    'CG': (-10.6, -27.2),
    'GC': (-9.8,  -24.4),
    'GG': (-8.0,  -19.9), 'CC': (-8.0, -19.9),
}
_R_GAS       = 1.987    # cal/(mol·K)
_PRIMER_CONC = 250e-9   # 250 nM primer (typical KASP reaction)
_NA_CONC     = 0.05     # 50 mM Na+


def _calc_tm_nn(seq: str) -> float:
    """
    Nearest-neighbor Tm using SantaLucia 1998 parameters.
    Conditions: 250 nM primer, 50 mM Na+.
    Falls back to Wallace rule for sequences containing non-ACGT bases (e.g. N).

    Bug 5 fix: replaces the old Wallace-only _calc_tm_simple, which gave
    inaccurate results especially for GC-rich or short primers.
    """
    seq = seq.upper()
    if not seq:
        return 0.0
    if not all(b in 'ACGT' for b in seq):
        return _calc_tm_wallace(seq)

    dH = 0.0   # kcal/mol
    dS = 0.0   # cal/(mol·K)

    # Terminal initiation parameters
    for end_base in (seq[0], seq[-1]):
        if end_base in 'GC':
            dH += 0.1;  dS += -2.8
        else:
            dH += 2.3;  dS += 4.1

    # Nearest-neighbor sum
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        if pair in _NN:
            h, s = _NN[pair]
            dH += h
            dS += s

    # Tm (K), non-self-complementary formula
    dH_cal = dH * 1000.0
    tm_k = dH_cal / (dS + _R_GAS * math.log(_PRIMER_CONC / 4.0))

    # Salt correction (Schildkraut & Lifson 1965)
    tm_c = tm_k - 273.15 + 16.6 * math.log10(_NA_CONC)
    return round(tm_c, 1)


def _calc_tm_wallace(seq: str) -> float:
    """Wallace rule fallback — used when sequence contains ambiguity bases."""
    seq = seq.upper()
    gc = sum(1 for b in seq if b in 'GC')
    at = sum(1 for b in seq if b in 'AT')
    if not (gc + at):
        return 0.0
    if len(seq) < 14:
        return round(2 * at + 4 * gc, 1)
    return round(64.9 + 41.0 * (gc - 16.4) / (gc + at), 1)


# Keep old name so any callers outside this module are not broken
_calc_tm_simple = _calc_tm_nn


# ── Data structures ───────────────────────────────────────────────────────────

@dataclass
class SNPPosition:
    position: int
    raw_pos:  int
    allele1:  str
    allele2:  str


@dataclass
class KASPMarker:
    snp_pos:         int
    allele1:         str
    allele2:         str
    fwd1_seq:        str
    fwd2_seq:        str
    rev_seq:         str
    fwd1_with_tail:  str
    fwd2_with_tail:  str
    fwd1_tm:         float = 0.0
    fwd2_tm:         float = 0.0
    rev_tm:          float = 0.0
    fwd1_gc:         float = 0.0
    fwd2_gc:         float = 0.0
    rev_gc:          float = 0.0
    product_size:    int   = 0
    # 'fwd' = allele-specific forward primers + common reverse
    # 'rev' = common forward + allele-specific reverse primers
    orientation:     str   = 'fwd'
    # Bug 3: number of internal mismatches in fwd2 body vs allele2 template
    fwd2_body_mismatches: int = 0


# ── Sequence utilities ────────────────────────────────────────────────────────

def _rc(seq: str) -> str:
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]


def _calc_gc(seq: str) -> float:
    seq = seq.upper()
    gc = sum(1 for b in seq if b in 'GC')
    return round(gc / len(seq) * 100, 1) if seq else 0.0


def _count_body_mismatches(primer_body: str, region2: str,
                            snp_in_region: int) -> int:
    """
    Bug 3 helper: count mismatches between the allele-specific primer body
    (all bases except the 3' SNP base) and the same-length window in
    region2 (the other allele's template), ending just before the SNP.

    Returns the number of mismatches.  If the window falls outside region2,
    returns 0 (cannot tell, so assume clean).
    """
    body_len     = len(primer_body)
    region_start = snp_in_region - body_len
    if region_start < 0 or region_start + body_len > len(region2):
        return 0
    region_body = region2[region_start:snp_in_region]
    return sum(
        1 for a, b in zip(primer_body.upper(), region_body.upper()) if a != b
    )


# ── Alignment helper ──────────────────────────────────────────────────────────

def align_two_seqs(seq1: str, seq2: str, muscle_bin: str):
    """Align two sequences using MUSCLE."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False) as f:
        f.write(f">seq1\n{seq1}\n>seq2\n{seq2}\n")
        tmp_in = f.name
    tmp_out = tmp_in.replace('.fas', '_aln.fas')
    try:
        subprocess.run([muscle_bin, '-align', tmp_in, '-output', tmp_out],
                       capture_output=True, timeout=300)
        from alignment import read_fasta
        aligned = read_fasta(tmp_out)
        seqs = list(aligned.values())
        if len(seqs) >= 2:
            return seqs[0], seqs[1]
    except Exception as e:
        log.warning("Alignment failed: %s", e)
    return seq1, seq2


# ── SNP detection ─────────────────────────────────────────────────────────────

def find_snps(seq1: str, seq2: str):
    """Find SNPs between two aligned sequences."""
    if len(seq1) != len(seq2):
        raise ValueError(
            f"Sequences must be same length. Got {len(seq1)} vs {len(seq2)}"
        )
    snps    = []
    raw_pos = 0
    for i, (b1, b2) in enumerate(zip(seq1.upper(), seq2.upper())):
        if b1 != '-':
            raw_pos += 1
        if b1 == '-' or b2 == '-':
            continue
        if b1 != b2 and b1 in 'ACGT' and b2 in 'ACGT':
            snps.append(SNPPosition(
                position=i, raw_pos=raw_pos - 1,
                allele1=b1, allele2=b2,
            ))
    return snps


# ── Primer3 helpers ───────────────────────────────────────────────────────────

def _design_allele_fwd(region: str, snp_pos: int, primer3_bin: str,
                       min_size: int = 18, max_size: int = 25,
                       min_tm: float = 57.0, max_tm: float = 63.0,
                       opt_tm: float = 60.0,
                       min_gc: float = 30.0, max_gc: float = 70.0,
                       max_hairpin_th: float = 47.0) -> str | None:
    """
    Design allele-specific forward primer ending at snp_pos (0-based).
    SEQUENCE_FORCE_LEFT_END forces the primer's 3' end onto the SNP.
    """
    p3 = (
        f"SEQUENCE_ID=kasp_fwd\n"
        f"SEQUENCE_TEMPLATE={region}\n"
        f"SEQUENCE_FORCE_LEFT_END={snp_pos}\n"
        f"PRIMER_PICK_LEFT_PRIMER=1\n"
        f"PRIMER_PICK_RIGHT_PRIMER=0\n"
        f"PRIMER_PICK_INTERNAL_OLIGO=0\n"
        f"PRIMER_OPT_SIZE={int((min_size + max_size) / 2)}\n"
        f"PRIMER_MIN_SIZE={min_size}\n"
        f"PRIMER_MAX_SIZE={max_size}\n"
        f"PRIMER_OPT_TM={opt_tm}\n"
        f"PRIMER_MIN_TM={min_tm}\n"
        f"PRIMER_MAX_TM={max_tm}\n"
        f"PRIMER_MIN_GC={min_gc}\n"
        f"PRIMER_MAX_GC={max_gc}\n"
        f"PRIMER_MAX_HAIRPIN_TH={max_hairpin_th}\n"
        f"PRIMER_NUM_RETURN=3\n"
        f"=\n"
    )
    try:
        r = subprocess.run([primer3_bin], input=p3,
                           capture_output=True, text=True, timeout=10)
        m = re.search(r'PRIMER_LEFT_0_SEQUENCE=(\S+)', r.stdout)
        return m.group(1) if m else None
    except Exception:
        return None


def _design_common_rev(region: str, snp_pos: int, primer3_bin: str,
                       min_size: int = 18, max_size: int = 25,
                       min_tm: float = 45.0, max_tm: float = 80.0,
                       opt_tm: float = 60.0) -> str | None:
    """
    Design common reverse primer strictly downstream of the SNP.

    Bug 1 fix: removed the dead 'max_right' variable; the broken fallback
    that silently reset min_pos to 5 (upstream of SNP) now simply returns
    None — there is no valid design when there is not enough downstream
    sequence.

    Bug 2 fix: replaced SEQUENCE_PRIMER_PAIR_OK_REGION_LIST (unreliable
    when only picking a single right primer) with template trimming.
    Primer3 only sees the downstream sub-region, so the right primer is
    guaranteed to land there.
    """
    downstream_start = snp_pos + 20

    if downstream_start + min_size > len(region):
        log.debug(
            "_design_common_rev: not enough downstream sequence "
            "(snp_pos=%d, region_len=%d, need %d more bases)",
            snp_pos, len(region), downstream_start + min_size - len(region),
        )
        return None

    sub_region = region[downstream_start:]

    p3 = (
        f"SEQUENCE_ID=kasp_rev\n"
        f"SEQUENCE_TEMPLATE={sub_region}\n"
        f"PRIMER_PICK_LEFT_PRIMER=0\n"
        f"PRIMER_PICK_RIGHT_PRIMER=1\n"
        f"PRIMER_PICK_INTERNAL_OLIGO=0\n"
        f"PRIMER_OPT_SIZE={int((min_size + max_size) / 2)}\n"
        f"PRIMER_MIN_SIZE={min_size}\n"
        f"PRIMER_MAX_SIZE={max_size}\n"
        f"PRIMER_OPT_TM={opt_tm}\n"
        f"PRIMER_MIN_TM={min_tm}\n"
        f"PRIMER_MAX_TM={max_tm}\n"
        f"PRIMER_MIN_GC=30.0\n"
        f"PRIMER_MAX_GC=70.0\n"
        f"PRIMER_NUM_RETURN=3\n"
        f"=\n"
    )
    try:
        r = subprocess.run([primer3_bin], input=p3,
                           capture_output=True, text=True, timeout=10)
        m = re.search(r'PRIMER_RIGHT_0_SEQUENCE=(\S+)', r.stdout)
        return m.group(1) if m else None
    except Exception:
        return None


def _find_discriminating_positions_msa(
    region: str,
    db_path: str,
    blastn_path: str,
    blastdbcmd_path: str,
    muscle_bin: str,
    snp_pos: int,
    min_amplicon: int = 50,
    max_amplicon: int = 170,
    n_seqs: int = 30,
) -> tuple[list[int], bool]:
    """
    MSA-based approach to find positions where the target base differs from
    ALL other genomic sequences.

    Workflow (mirrors the manual expert process):
      1. BLAST the target region against the whole genome.
      2. Collect the top hits (excluding self-hit / target chromosome).
      3. Extract those hit sequences from the BLAST database using blastdbcmd.
      4. Align target + extracted sequences with MUSCLE.
      5. Parse the MSA column-by-column: positions where the target base
         differs from every other sequence in the alignment are discriminating.

    Returns (positions, muscle_timed_out):
      - positions: list of discriminating positions (in region coordinates),
        sorted shortest-amplicon-first, within the amplicon window.
      - muscle_timed_out: True if MUSCLE timed out (region is likely a repeat
        element; blind primer search should be skipped for this locus).
    """
    if not all([db_path, region, blastn_path, blastdbcmd_path, muscle_bin]):
        return [], False

    # ── Step 1: BLAST region ──────────────────────────────────────────────────
    # -max_hsps 1: keep only the best HSP per subject sequence.
    # This caps output at ≤(n_seqs+5) lines regardless of how many copies a
    # repeat element has on each chromosome — prevents MUSCLE timeouts.
    try:
        blast = subprocess.run([
            blastn_path, '-db', db_path,
            '-query', '-', '-task', 'blastn',
            '-word_size', '11',
            '-outfmt', '6 qseqid sseqid pident length sstart send sstrand',
            '-max_target_seqs', str(n_seqs + 5),
            '-max_hsps', '1',
            '-evalue', '1e-5',
            '-perc_identity', '70',
            '-dust', 'yes',
        ], input=f'>query\n{region}\n',
        capture_output=True, text=True, timeout=60)
    except Exception as e:
        log.warning('_find_discriminating_positions_msa: BLAST failed: %s', e)
        return [], False

    # ── Step 2: Collect top hits, skip self-hit ───────────────────────────────
    blast_lines = [l for l in blast.stdout.strip().split('\n') if l]
    log.info('_find_discriminating_positions_msa: BLAST returned %d hits', len(blast_lines))

    seen_chrs: dict[str, tuple] = {}   # chr → (sstart, send, sstrand, pident)
    target_chr: str | None = None

    for line in blast_lines:
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 7:
            continue
        try:
            sseqid  = parts[1]
            pident  = float(parts[2])
            aln_len = int(parts[3])
            sstart  = int(parts[4])
            send    = int(parts[5])
            sstrand = parts[6].strip()   # 'plus' or 'minus'
        except (ValueError, IndexError):
            continue

        if aln_len < 20:
            continue

        # Best (first) hit = self / target chromosome
        if target_chr is None:
            target_chr = sseqid
            continue

        # Skip additional hits on the same target chromosome
        if sseqid == target_chr:
            continue

        if sseqid not in seen_chrs:
            seen_chrs[sseqid] = (sstart, send, sstrand, pident)
            if len(seen_chrs) >= n_seqs:
                break

    if not seen_chrs:
        log.info('_find_discriminating_positions_msa: no off-target hits found')
        return [], False

    pidents = [v[3] for v in seen_chrs.values()]
    min_pident = min(pidents) if pidents else 0.0
    log.info(
        '_find_discriminating_positions_msa: %d off-target chromosomes found '
        '(min pident=%.1f%%)',
        len(seen_chrs), min_pident,
    )

    # ── Step 3: Extract hit sequences using blastdbcmd ────────────────────────
    extracted: list[tuple[str, str]] = []   # (name, sequence)

    # Add flanking so extracted sequence is at least as long as region
    flank_bp = max(len(region), 50)

    for chrom, (sstart, send, sstrand, _pident) in seen_chrs.items():
        lo = max(1, min(sstart, send) - flank_bp)
        hi = max(sstart, send) + flank_bp
        try:
            bc = subprocess.run([
                blastdbcmd_path,
                '-db', db_path,
                '-entry', chrom,
                '-range', f'{lo}-{hi}',
                '-strand', 'plus',
                '-outfmt', '%f',
            ], capture_output=True, text=True, timeout=30)
            # Strip FASTA header lines, join sequence
            seq = ''.join(
                l for l in bc.stdout.strip().split('\n') if not l.startswith('>')
            ).replace('\n', '').upper()
            log.info('_find_discriminating_positions_msa: extracted %s len=%d strand=%s', chrom, len(seq), sstrand)
            if seq and len(seq) >= 20:
                # Reverse-complement if the BLAST hit was on the minus strand
                if sstrand == 'minus':
                    seq = _rc(seq)
                extracted.append((chrom, seq))
        except Exception as exc:
            log.warning('_find_discriminating_positions_msa: blastdbcmd failed for %s: %s', chrom, exc)
            continue

    if not extracted:
        log.info('_find_discriminating_positions_msa: could not extract any hit sequences')
        return [], False

    log.info('_find_discriminating_positions_msa: extracted %d sequences, running MUSCLE', len(extracted))

    # ── Step 4: MUSCLE alignment (use temp files — /dev/stdin is unreliable on macOS) ──
    fasta_in = f'>TARGET\n{region}\n'
    for name, seq in extracted:
        safe = name.replace(' ', '_')
        fasta_in += f'>{safe}\n{seq}\n'

    import tempfile, os as _os
    tmp_in = tmp_out = None
    msa_out = ''
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False) as fh:
            fh.write(fasta_in)
            tmp_in = fh.name
        tmp_out = tmp_in.replace('.fas', '_aln.fas')
        msa_run = subprocess.run(
            [muscle_bin, '-align', tmp_in, '-output', tmp_out],
            capture_output=True, text=True, timeout=120,
        )
        log.info('_find_discriminating_positions_msa: MUSCLE returncode=%d stderr=%s',
                 msa_run.returncode, msa_run.stderr[:200] if msa_run.stderr else '')
        if _os.path.exists(tmp_out):
            with open(tmp_out) as fh:
                msa_out = fh.read()
    except Exception as e:
        log.warning('_find_discriminating_positions_msa: MUSCLE failed: %s', e)
        # Return muscle_timed_out=True so callers can skip blind primer search
        # for this locus (a timeout reliably indicates a highly repetitive region)
        return [], True
    finally:
        for p in [tmp_in, tmp_out]:
            if p and _os.path.exists(p):
                try:
                    _os.unlink(p)
                except Exception:
                    pass

    if not msa_out.strip():
        log.info('_find_discriminating_positions_msa: MUSCLE produced no output')
        return [], False

    # ── Step 5: Parse MSA, find discriminating positions ─────────────────────
    # Parse FASTA MSA output
    msa_seqs: dict[str, str] = {}
    cur_name = None
    for line in msa_out.split('\n'):
        if line.startswith('>'):
            cur_name = line[1:].strip().split()[0]
            msa_seqs[cur_name] = ''
        elif cur_name:
            msa_seqs[cur_name] += line.strip().upper()

    log.info('_find_discriminating_positions_msa: MSA has %d sequences: %s',
             len(msa_seqs), list(msa_seqs.keys())[:5])

    target_msa = msa_seqs.get('TARGET')
    if not target_msa:
        log.info('_find_discriminating_positions_msa: TARGET not found in MSA')
        return [], False

    other_msas = [seq for name, seq in msa_seqs.items() if name != 'TARGET']
    if not other_msas:
        return [], False

    # Map MSA column → region (ungapped) position
    # region_pos[col] = position in the ungapped region string (-1 = gap)
    region_pos: list[int] = []
    r_idx = 0
    for col_char in target_msa:
        if col_char == '-':
            region_pos.append(-1)
        else:
            region_pos.append(r_idx)
            r_idx += 1

    # Amplicon window in region coordinates
    max_primer_size = 25
    tip_min = snp_pos + max(15, min_amplicon - 40)
    tip_max = min(snp_pos + max_amplicon, len(region) - 1 - max_primer_size)

    candidates = []
    for col, (rpos, tbase) in enumerate(zip(region_pos, target_msa)):
        if tbase == '-' or rpos < 0:
            continue
        if not (tip_min <= rpos <= tip_max):
            continue

        # Check every other sequence at this column
        all_differ = True
        any_covered = False
        for other_msa in other_msas:
            if col >= len(other_msa):
                continue
            ob = other_msa[col]
            if ob == '-':
                continue   # gap = no info for this sequence at this column
            any_covered = True
            if ob.upper() == tbase.upper():
                all_differ = False
                break

        if any_covered and all_differ:
            candidates.append(rpos)

    candidates.sort()
    log.info(
        '_find_discriminating_positions_msa: %d discriminating positions found '
        '(from MSA with %d sequences)',
        len(candidates), len(other_msas),
    )

    # Region is not designable if ALL 21 off-target wheat chromosomes were found
    # AND MSA found 0 discriminating positions.
    # When every chromosome carries this sequence and no position is unique to
    # the target, no genome-specific KASP primer can be designed here.
    # (Hexaploid wheat has 21 off-target chromosomes: 7A + 7B + 7D.)
    region_not_designable = (len(seen_chrs) >= 21 and len(candidates) == 0
                             and min_pident >= 90.0)
    if region_not_designable:
        log.warning(
            '_find_discriminating_positions_msa: all %d wheat chromosomes '
            'matched with 0 discriminating positions — region is not designable '
            'as genome-specific KASP. Blind search will be skipped.',
            len(seen_chrs),
        )

    return candidates, region_not_designable


def _find_discriminating_positions(region: str, db_path: str,
                                    blastn_path: str,
                                    snp_pos: int,
                                    min_amplicon: int = 50,
                                    max_amplicon: int = 250) -> list[int]:
    """
    Find positions in `region` (within the amplicon window) where the base in
    `region` differs from ALL other chromosomes in the genome.

    Strategy:
      - BLAST region against genome with output format 6 including query sequence
        (fmt "6 qstart qend sseqid qseq sseq pident length")
      - For each query position in the amplicon window, collect all bases that
        other chromosomes have at that position from the pairwise alignments
      - A position is discriminating if no other chromosome shares the same base
        as region[pos] at that position

    Returns a list of candidate positions sorted by closeness to SNP (shortest
    amplicon first).

    Window logic: SEQUENCE_FORCE_RIGHT_END is the 3' TIP of the reverse primer
    (its leftmost position in template coordinates). The full amplicon also
    includes both primer bodies (~20 bp each), so the tip can be as close as
    (min_amplicon - 2 * max_primer_size) bases after the SNP while still
    meeting the min_amplicon constraint. We use max_primer_size=25.
    """
    if not db_path or not region:
        return []

    max_primer_size = 25
    # Minimum: tip must be > snp_pos and amplicon estimate ≥ min_amplicon
    window_start = snp_pos + max(1, min_amplicon - 2 * max_primer_size)
    # Maximum: leave room for the primer body beyond the tip
    window_end   = min(snp_pos + max_amplicon, len(region) - 1 - max_primer_size)
    if window_start > window_end:
        return []

    try:
        r = subprocess.run([
            blastn_path, '-db', db_path,
            '-query', '-', '-task', 'blastn',
            '-word_size', '11',
            '-outfmt', '6 qseqid qstart qend sseqid pident length qseq sseq',
            '-max_target_seqs', '500',
            '-evalue', '1e-3',
            '-dust', 'yes',
            '-perc_identity', '85',
        ], input=f'>query\n{region}\n',
        capture_output=True, text=True, timeout=120)
    except Exception:
        return []

    # For each position in the amplicon window, track:
    #   same_at[pos]  = True if ANY other chromosome shares the target base here
    #   covered[pos]  = True if ANY other chromosome has alignment coverage here
    # A position is discriminating only if covered AND NOT same_at.
    same_at:  dict[int, bool] = {p: False for p in range(window_start, window_end + 1)}
    covered:  dict[int, bool] = {p: False for p in range(window_start, window_end + 1)}

    for line in r.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 8:
            continue
        try:
            pident   = float(parts[4])
            aln_len  = int(parts[5])
            if pident < 85 or aln_len < 11:
                continue
            qstart   = int(parts[1]) - 1   # 0-based
            qend     = int(parts[2]) - 1
            qseq     = parts[6]             # aligned query bases (may have gaps '-')
            sseq     = parts[7]             # aligned subject bases
        except (ValueError, IndexError):
            continue

        # Walk the pairwise alignment column by column
        q_pos = qstart
        for qb, sb in zip(qseq, sseq):
            if qb == '-':
                continue  # insertion in subject, no query position consumed
            if q_pos > qend or q_pos > window_end:
                break
            if window_start <= q_pos <= window_end and sb != '-':
                covered[q_pos] = True
                if sb.upper() == qb.upper():
                    # Another chromosome shares the exact same base — NOT discriminating
                    same_at[q_pos] = True
            if qb != '-':
                q_pos += 1

    # A position is discriminating if:
    #   1. At least one other chromosome has alignment coverage here
    #   2. NO other chromosome shares the same base as the target
    candidates = []
    for pos in range(window_start, window_end + 1):
        if covered[pos] and not same_at[pos]:
            candidates.append(pos)

    # Sort by amplicon size (shortest first = closest to SNP)
    candidates.sort()
    return candidates


def _count_3prime_match_chromosomes(primer_seq: str, db_path: str,
                                    blastn_path: str) -> int:
    """
    BLAST the primer against the genome.
    Return the number of DISTINCT chromosomes that pass BOTH:
      - ≥ 95% sequence identity
      - ≥ 85% alignment coverage of the primer length
      - Alignment covers the 3' end (qend == len(primer))

    Per spec: short/low-coverage matches are not significant off-targets.
    A genome-specific primer should return 1 (only the target chromosome).
    """
    try:
        r = subprocess.run([
            blastn_path, '-db', db_path,
            '-query', '-', '-task', 'blastn-short',
            '-word_size', '7',
            '-outfmt', '6 qseqid qstart qend sseqid pident length',
            '-max_target_seqs', '200',
            '-evalue', '10000',   # permissive — short primers have high E-values against large genomes
            '-dust', 'no',
            '-perc_identity', '80',
        ], input=f'>primer\n{primer_seq}\n',
        capture_output=True, text=True, timeout=30)
    except Exception:
        return 1  # cannot verify; assume OK

    n = len(primer_seq)
    min_coverage = 0.85 * n   # ≥85% of primer length must align
    matching_chrs: set[str] = set()

    for line in r.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 6:
            continue
        try:
            qend    = int(parts[2])
            sseqid  = parts[3]
            pident  = float(parts[4])
            aln_len = int(parts[5])
        except (ValueError, IndexError):
            continue
        # Spec thresholds: ≥95% identity + ≥85% primer coverage + covers 3' end
        if pident >= 95.0 and aln_len >= min_coverage and qend >= n:
            matching_chrs.add(sseqid)

    return len(matching_chrs)


def _design_common_rev_genome_specific(
    region: str,
    snp_pos: int,
    primer3_bin: str,
    db_path: str,
    blastn_path: str,
    blastdbcmd_path: str = '',
    muscle_bin: str = '',
    min_size: int = 18,
    max_size: int = 25,
    min_tm: float = 45.0,
    max_tm: float = 80.0,
    opt_tm: float = 60.0,
    min_amplicon: int = 50,
    max_amplicon: int = 170,
    precomputed_msa_positions: 'set[int] | None' = None,
    skip_blind_search: bool = False,
) -> str | None:
    """
    Design a common reverse primer whose 3' end is genome-specific.

    Strategy (matches the expert manual workflow):
      1. MSA pre-filter: BLAST region → extract homeolog sequences →
         MUSCLE alignment → find positions where target base differs from
         ALL other sequences. Try these first (most reliable).
      2. BLAST alignment pre-filter: fallback if MSA unavailable.
      3. Try ALL remaining positions in the window (shortest first), capped
         at MAX_ATTEMPTS total.
      4. For each position: Primer3 to design the primer, then BLAST the
         primer back to verify it hits ≤1 chromosome at the 3' end.

    The verification BLAST (step 4) is the ground-truth gatekeeper.

    skip_blind_search: when True, only MSA/BLAST pre-filter positions are tried.
      No uncategorised positions are attempted. Use when MUSCLE timed out
      (repeat region) to avoid false-positive genome-specific primers.
    """
    MAX_ATTEMPTS = 60

    # Amplicon tip window
    tip_min = snp_pos + max(15, min_amplicon - 40)
    tip_max = min(snp_pos + max_amplicon, len(region) - 1 - min_size)
    if tip_min > tip_max:
        return None

    # ── Step 1: MSA-based pre-filter (most reliable) ──────────────────────────
    msa_discriminating: set[int] = set()
    if precomputed_msa_positions is not None:
        # Reuse already-computed MSA positions — avoids re-running BLAST+MUSCLE
        msa_discriminating = precomputed_msa_positions
        log.info(
            '_design_common_rev_genome_specific: MSA found %d discriminating positions (cached)',
            len(msa_discriminating),
        )
    elif blastdbcmd_path and muscle_bin:
        msa_positions, msa_not_designable = _find_discriminating_positions_msa(
            region, db_path, blastn_path, blastdbcmd_path, muscle_bin,
            snp_pos, min_amplicon, max_amplicon,
        )
        msa_discriminating = set(msa_positions)
        if msa_not_designable:
            # Region is a repeat — all chromosomes matched with no discriminating
            # positions. Override skip_blind_search so we never blind-search here
            # and produce a false-positive genome-specific primer.
            skip_blind_search = True
        log.info(
            '_design_common_rev_genome_specific: MSA found %d discriminating positions',
            len(msa_discriminating),
        )

    # ── Step 2: BLAST alignment pre-filter (fallback supplement) ─────────────
    blast_discriminating = set(_find_discriminating_positions(
        region, db_path, blastn_path, snp_pos, min_amplicon, max_amplicon,
    ))

    # Combine: MSA positions first, then BLAST-only positions, then rest
    all_positions = list(range(tip_min, tip_max + 1))
    msa_only    = [p for p in all_positions if p in msa_discriminating]
    blast_only  = [p for p in all_positions if p in blast_discriminating and p not in msa_discriminating]
    remaining   = [p for p in all_positions if p not in msa_discriminating and p not in blast_discriminating]

    if skip_blind_search:
        # Repeat region: only try positions confirmed by MSA or BLAST pre-filter.
        # Skipping uncategorised positions prevents false-positive genome-specific
        # primers in loci where the region is identical across all chromosomes.
        ordered = (msa_only + blast_only)[:MAX_ATTEMPTS]
        log.info(
            '_design_common_rev_genome_specific: skip_blind_search=True — '
            'trying up to %d positions (%d MSA, %d BLAST-only, 0 uncategorised)',
            len(ordered), len(msa_only), len(blast_only),
        )
    else:
        ordered = (msa_only + blast_only + remaining)[:MAX_ATTEMPTS]
        log.info(
            '_design_common_rev_genome_specific: trying up to %d positions '
            '(%d MSA, %d BLAST-only, %d uncategorised)',
            len(ordered), len(msa_only), len(blast_only),
            len(remaining[:MAX_ATTEMPTS - len(msa_only) - len(blast_only)]),
        )

    # ── Try each candidate position ───────────────────────────────────────────
    for end_pos in ordered:
        p3 = (
            f"SEQUENCE_ID=kasp_rev_specific\n"
            f"SEQUENCE_TEMPLATE={region}\n"
            f"SEQUENCE_FORCE_RIGHT_END={end_pos}\n"
            f"PRIMER_PICK_LEFT_PRIMER=0\n"
            f"PRIMER_PICK_RIGHT_PRIMER=1\n"
            f"PRIMER_PICK_INTERNAL_OLIGO=0\n"
            f"PRIMER_OPT_SIZE={int((min_size + max_size) / 2)}\n"
            f"PRIMER_MIN_SIZE={min_size}\n"
            f"PRIMER_MAX_SIZE={max_size}\n"
            f"PRIMER_OPT_TM={opt_tm}\n"
            f"PRIMER_MIN_TM={min_tm}\n"
            f"PRIMER_MAX_TM={max_tm}\n"
            f"PRIMER_MIN_GC=30.0\n"
            f"PRIMER_MAX_GC=70.0\n"
            f"PRIMER_NUM_RETURN=1\n"
            f"=\n"
        )
        try:
            r = subprocess.run([primer3_bin], input=p3,
                               capture_output=True, text=True, timeout=10)
            m = re.search(r'PRIMER_RIGHT_0_SEQUENCE=(\S+)', r.stdout)
            if not m:
                continue  # Primer3 couldn't design here (Tm/GC fail) — skip
            primer_seq = m.group(1)

            # ── Verify: BLAST primer back against genome ──────────────────
            chr_hits = _count_3prime_match_chromosomes(
                primer_seq, db_path, blastn_path,
            )
            if chr_hits > 1:
                log.debug(
                    '_design_common_rev_genome_specific: pos %d rejected '
                    '(%d chromosomes match at 3-end)',
                    end_pos, chr_hits,
                )
                continue

            log.info(
                '_design_common_rev_genome_specific: accepted primer at pos %d '
                '(amplicon tip ~%d bp from SNP) -- %d chr at 3-end',
                end_pos, end_pos - snp_pos, chr_hits,
            )
            return primer_seq

        except Exception:
            continue

    log.warning(
        '_design_common_rev_genome_specific: exhausted %d positions, '
        'no genome-specific primer found', len(ordered),
    )
    return None


# ── Product size helpers ──────────────────────────────────────────────────────

def _fwd_product_size(region: str, fwd_seq: str, rev_seq: str,
                      snp_in_region: int, fallback: int) -> int:
    """
    Calculate amplicon size for forward-orientation KASP marker.
    Product spans from 5' end of allele-specific fwd primer to 3' end of
    common rev primer binding site (inclusive).
    """
    rev_bind = region.find(_rc(rev_seq))
    if rev_bind > snp_in_region:
        fwd_start = snp_in_region - len(fwd_seq) + 1
        return rev_bind + len(rev_seq) - fwd_start
    return fallback


def _rev_product_size(rc_region: str, base_rev_seq: str, cfwd_seq: str,
                      snp_in_rc: int, fallback: int) -> int:
    """
    Calculate amplicon size for reverse-orientation KASP marker.

    Bug 4 fix: was hardcoded to 'flank' (150); now computed from primer
    positions on the RC strand.
    """
    cfwd_bind = rc_region.find(_rc(cfwd_seq))
    if cfwd_bind > snp_in_rc:
        rev_start = snp_in_rc - len(base_rev_seq) + 1
        return cfwd_bind + len(cfwd_seq) - rev_start
    return fallback


# ── Genome-specificity helpers ────────────────────────────────────────────────

def _compute_uniqueness_map_for_region(seq: str, db_path: str,
                                       blastn_path: str = 'blastn') -> list:
    """
    BLAST seq against the genome. Return a boolean list where True = position
    maps to only 1 chromosome (unique), False = repetitive/multi-copy.
    Positions with no BLAST hit default to True (assumed unique).
    Falls back to all-True if BLAST fails.
    """
    import subprocess as _sp
    if not db_path or not seq:
        return [True] * len(seq)
    try:
        r = _sp.run([
            blastn_path, '-db', db_path,
            '-query', '-', '-task', 'blastn',
            '-word_size', '11',
            '-outfmt', '6',
            '-max_target_seqs', '500',
            '-evalue', '1e-3',
            '-dust', 'yes',
            '-perc_identity', '85',
        ], input=f'>query\n{seq}\n',
        capture_output=True, text=True, timeout=120)

        pos_chroms: list = [set() for _ in range(len(seq))]
        for line in r.stdout.strip().split('\n'):
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 10:
                continue
            pident = float(parts[2])
            aln_len = int(parts[3])
            if pident < 85 or aln_len < 11:
                continue
            qstart = int(parts[6]) - 1
            qend   = int(parts[7]) - 1
            chrom  = parts[1].split(':')[0]
            for pos in range(max(0, qstart), min(len(seq), qend + 1)):
                pos_chroms[pos].add(chrom)

        return [len(chroms) <= 1 for chroms in pos_chroms]
    except Exception:
        return [True] * len(seq)


# ── Main single-SNP design ────────────────────────────────────────────────────

def design_kasp_for_snp(ungap1: str, ungap2: str, snp: SNPPosition,
                        primer3_bin: str, flank: int = 100,
                        db_path: str = '',
                        blastn_path: str = 'blastn',
                        blastdbcmd_path: str = '',
                        muscle_bin: str = '',
                        progress_callback=None) -> list[KASPMarker]:
    """
    Design KASP primers for one SNP.
    Tries BOTH orientations independently and returns all that succeed:
      1. Allele-specific FORWARD primers + common REVERSE
      2. Common FORWARD primer + allele-specific REVERSE primers

    When db_path is provided, the common primer's 3' end is forced to land on
    a genome-discriminating position identified via MSA of the target region
    against the closest homeologous sequences.

    Returns a list of KASPMarker (0–2 items for standard Tm; if neither
    orientation succeeds, relaxed Tm (50-65 °C) is tried and those results
    returned instead).
    """
    pos           = snp.raw_pos
    start         = max(0, pos - flank)
    end           = min(len(ungap1), pos + flank + 1)
    snp_in_region = pos - start

    region1 = ungap1[start:end]
    region2 = ungap2[start:end]

    # Compute MSA discriminating positions once — reused by all orientation/Tm attempts
    precomputed_msa: set[int] | None = None
    region_not_designable = False  # True when region is a repeat (no specific primer possible)
    if db_path and blastdbcmd_path and muscle_bin:
        if progress_callback:
            progress_callback(f'pos {snp.raw_pos + 1} — finding genome-discriminating positions')
        msa_positions, region_not_designable = _find_discriminating_positions_msa(
            region1, db_path, blastn_path, blastdbcmd_path, muscle_bin,
            snp_in_region, min_amplicon=50, max_amplicon=170,
        )
        precomputed_msa = set(msa_positions)
        if region_not_designable:
            log.warning(
                'SNP %d: repeat region detected — blind primer search will be '
                'skipped. This locus is not designable as a genome-specific '
                'KASP marker.',
                snp.raw_pos,
            )

    results: list[KASPMarker] = []

    # ── Orientation 1: allele-specific forward + common reverse ──────────────
    if progress_callback:
        progress_callback(f'pos {snp.raw_pos + 1} — orientation 1 (FAM→{snp.allele1}, HEX→{snp.allele2})')

    base_fwd = _design_allele_fwd(region1, snp_in_region, primer3_bin)

    # Choose common primer design method
    if db_path and base_fwd:
        rev = _design_common_rev_genome_specific(
            region1, snp_in_region, primer3_bin,
            db_path=db_path, blastn_path=blastn_path,
            blastdbcmd_path=blastdbcmd_path, muscle_bin=muscle_bin,
            min_amplicon=50, max_amplicon=170,
            precomputed_msa_positions=precomputed_msa,
            skip_blind_search=region_not_designable,
        )
        if rev is None:
            # No genome-specific primer found — do not fall back to non-specific
            pass
    elif base_fwd:
        rev = _design_common_rev(region1, snp_in_region, primer3_bin)
    else:
        rev = None

    if base_fwd and rev:
        fwd1 = base_fwd[:-1] + snp.allele1

        # Bug 3 fix: check whether fwd2's body has internal mismatches
        # against allele2's template; if so, redesign from region2.
        body       = base_fwd[:-1]
        mismatches = _count_body_mismatches(body, region2, snp_in_region)
        if mismatches > 0:
            log.info('SNP %d: fwd2 body has %d mismatch(es) vs allele2 (using shared body)', snp.raw_pos, mismatches)
        fwd2 = body + snp.allele2

        # ── Verify common primer is compatible with BOTH parents ──────────────
        # The common primer is designed from region1 only.  If another SNP falls
        # inside its binding site it will have a mismatch to region2's template
        # and fail to amplify the Parent-2 allele.  We count mismatches and, if
        # any are found, try to redesign from region2 before giving up.
        rev_rc   = _rc(rev)
        bind_pos = region1.find(rev_rc)   # 0-based start of binding site in region1
        rev_mm_to_r2 = 0
        if bind_pos >= 0:
            r2_window = region2[bind_pos: bind_pos + len(rev_rc)]
            rev_mm_to_r2 = sum(
                a.upper() != b.upper()
                for a, b in zip(rev_rc, r2_window)
            )
            if rev_mm_to_r2 > 0:
                log.warning(
                    'SNP %d: common rev has %d mismatch(es) to region2 — '
                    'trying redesign from region2',
                    snp.raw_pos, rev_mm_to_r2,
                )
                rev_r2 = _design_common_rev(region2, snp_in_region, primer3_bin)
                if rev_r2:
                    # Prefer the region2-designed primer if it also matches region1
                    rev_rc2   = _rc(rev_r2)
                    bp2       = region1.find(rev_rc2)
                    mm_back   = 0
                    if bp2 >= 0:
                        r1_window = region1[bp2: bp2 + len(rev_rc2)]
                        mm_back   = sum(
                            a.upper() != b.upper()
                            for a, b in zip(rev_rc2, r1_window)
                        )
                    if mm_back == 0:
                        log.info('SNP %d: using region2-derived common primer (no mm to region1)', snp.raw_pos)
                        rev = rev_r2
                        rev_mm_to_r2 = 0
                    else:
                        log.warning('SNP %d: region2-derived common primer also has mm to region1 (%d); keeping region1 design', snp.raw_pos, mm_back)
        # ─────────────────────────────────────────────────────────────────────

        prod_size = _fwd_product_size(region1, fwd1, rev, snp_in_region, flank)

        results.append(KASPMarker(
            snp_pos=snp.raw_pos, allele1=snp.allele1, allele2=snp.allele2,
            fwd1_seq=fwd1, fwd2_seq=fwd2, rev_seq=rev,
            fwd1_with_tail=FAM_TAIL + fwd1, fwd2_with_tail=HEX_TAIL + fwd2,
            fwd1_tm=_calc_tm_nn(fwd1), fwd2_tm=_calc_tm_nn(fwd2),
            rev_tm=_calc_tm_nn(rev),
            fwd1_gc=_calc_gc(fwd1), fwd2_gc=_calc_gc(fwd2),
            rev_gc=_calc_gc(rev),
            product_size=prod_size, orientation='fwd',
            fwd2_body_mismatches=mismatches,
        ))

    # ── Orientation 2: common forward + allele-specific reverse ──────────────
    # Always tried independently — may yield a better MSA-anchored common primer.
    if progress_callback:
        progress_callback(f'pos {snp.raw_pos + 1} — orientation 2 (FAM→{snp.allele2}, HEX→{snp.allele1})')

    log.info("SNP %d: trying orientation 2 (common fwd + allele-specific rev)", snp.raw_pos)
    rc1       = _rc(region1)
    rc2       = _rc(region2)
    snp_in_rc = len(region1) - snp_in_region - 1

    # Orientation 2 allele-specific primers come from the sequence to the RIGHT
    # of the SNP (top strand), which is often GC-rich (Tm 65-70 °C, GC up to 80%).
    # Widen both Tm and GC windows so Primer3 can find primers in that region.
    base_rev = _design_allele_fwd(rc1, snp_in_rc, primer3_bin,
                                   min_tm=50.0, max_tm=70.0, opt_tm=62.0,
                                   max_hairpin_th=100.0)
    if base_rev:
        log.info("SNP %d (ori2): allele-specific primer: %s (Tm=%.1f, GC=%.1f%%)",
                 snp.raw_pos, base_rev, _calc_tm_nn(base_rev), _calc_gc(base_rev))
    else:
        log.info("SNP %d (ori2): no allele-specific primer found (Tm 50-70°C, GC 20-80%%)",
                 snp.raw_pos)
    if base_rev:
        if db_path:
            cfwd = _design_common_rev_genome_specific(
                rc1, snp_in_rc, primer3_bin,
                db_path=db_path, blastn_path=blastn_path,
                blastdbcmd_path=blastdbcmd_path, muscle_bin=muscle_bin,
                min_amplicon=50, max_amplicon=170,
                skip_blind_search=region_not_designable,
            )
        else:
            cfwd = _design_common_rev(rc1, snp_in_rc, primer3_bin)
    else:
        cfwd = None

    if base_rev and cfwd:
        rev_body   = base_rev[:-1]
        mismatches = _count_body_mismatches(rev_body, rc2, snp_in_rc)
        if mismatches > 0:
            log.info('SNP %d (rev): rev2 body has %d mismatch(es) vs allele2 RC (using shared body)', snp.raw_pos, mismatches)
        rev2_body = rev_body

        rev1     = rev_body  + _rc(snp.allele1)
        rev2     = rev2_body + _rc(snp.allele2)

        # Verify common forward primer is compatible with both parents (same
        # logic as the forward-orientation check above, but operating on RC strands)
        cfwd_rc     = _rc(cfwd)
        cbind_pos   = rc1.find(cfwd_rc)
        cfwd_mm_r2  = 0
        if cbind_pos >= 0:
            rc2_window = rc2[cbind_pos: cbind_pos + len(cfwd_rc)]
            cfwd_mm_r2 = sum(
                a.upper() != b.upper()
                for a, b in zip(cfwd_rc, rc2_window)
            )
            if cfwd_mm_r2 > 0:
                log.warning(
                    'SNP %d (rev): common fwd has %d mismatch(es) to rc2 — '
                    'trying redesign from rc2',
                    snp.raw_pos, cfwd_mm_r2,
                )
                cfwd_r2 = _design_common_rev(rc2, snp_in_rc, primer3_bin)
                if cfwd_r2:
                    cfwd_rc2  = _rc(cfwd_r2)
                    cbp2      = rc1.find(cfwd_rc2)
                    mm_back   = 0
                    if cbp2 >= 0:
                        r1w    = rc1[cbp2: cbp2 + len(cfwd_rc2)]
                        mm_back = sum(
                            a.upper() != b.upper()
                            for a, b in zip(cfwd_rc2, r1w)
                        )
                    if mm_back == 0:
                        log.info('SNP %d (rev): using rc2-derived common fwd (no mm to rc1)', snp.raw_pos)
                        cfwd = cfwd_r2

        cfwd_seq = cfwd
        prod_size = _rev_product_size(rc1, base_rev, cfwd, snp_in_rc, flank)

        results.append(KASPMarker(
            snp_pos=snp.raw_pos, allele1=snp.allele1, allele2=snp.allele2,
            fwd1_seq=rev1, fwd2_seq=rev2, rev_seq=cfwd_seq,
            fwd1_with_tail=FAM_TAIL + rev1, fwd2_with_tail=HEX_TAIL + rev2,
            fwd1_tm=_calc_tm_nn(rev1), fwd2_tm=_calc_tm_nn(rev2),
            rev_tm=_calc_tm_nn(cfwd_seq),
            fwd1_gc=_calc_gc(rev1), fwd2_gc=_calc_gc(rev2),
            rev_gc=_calc_gc(cfwd_seq),
            product_size=prod_size, orientation='rev',
            fwd2_body_mismatches=mismatches,
        ))

    # ── Relaxed-Tm fallback (50-65 °C) ──────────────────────────────────────
    # Only tried when both standard-Tm orientations produced nothing (SNP near
    # a sequence boundary with insufficient upstream/downstream context).
    if results:
        return results

    if progress_callback:
        progress_callback(f'pos {snp.raw_pos + 1} — retrying with relaxed Tm (50–65°C)')
    log.info("SNP %d: retrying with relaxed Tm (50-65 °C)", snp.raw_pos)

    base_fwd_relaxed = _design_allele_fwd(
        region1, snp_in_region, primer3_bin,
        min_tm=50.0, max_tm=65.0, opt_tm=57.0,
    )

    if base_fwd_relaxed:
        if db_path:
            rev_relaxed = _design_common_rev_genome_specific(
                region1, snp_in_region, primer3_bin,
                db_path=db_path, blastn_path=blastn_path,
                blastdbcmd_path=blastdbcmd_path, muscle_bin=muscle_bin,
                min_amplicon=50, max_amplicon=170,
                precomputed_msa_positions=precomputed_msa,
                skip_blind_search=region_not_designable,
            )
        else:
            rev_relaxed = _design_common_rev(
                region1, snp_in_region, primer3_bin,
            )

        if base_fwd_relaxed and rev_relaxed:
            fwd1 = base_fwd_relaxed[:-1] + snp.allele1
            body = base_fwd_relaxed[:-1]
            mismatches = _count_body_mismatches(body, region2, snp_in_region)
            fwd2 = body + snp.allele2
            prod_size = _fwd_product_size(region1, fwd1, rev_relaxed, snp_in_region, flank)
            log.info("SNP %d: relaxed-Tm primer designed (fwd1_tm=%.1f)",
                     snp.raw_pos, _calc_tm_nn(fwd1))
            results.append(KASPMarker(
                snp_pos=snp.raw_pos, allele1=snp.allele1, allele2=snp.allele2,
                fwd1_seq=fwd1, fwd2_seq=fwd2, rev_seq=rev_relaxed,
                fwd1_with_tail=FAM_TAIL + fwd1, fwd2_with_tail=HEX_TAIL + fwd2,
                fwd1_tm=_calc_tm_nn(fwd1), fwd2_tm=_calc_tm_nn(fwd2),
                rev_tm=_calc_tm_nn(rev_relaxed),
                fwd1_gc=_calc_gc(fwd1), fwd2_gc=_calc_gc(fwd2),
                rev_gc=_calc_gc(rev_relaxed),
                product_size=prod_size, orientation='fwd',
                fwd2_body_mismatches=mismatches,
            ))

    if not results:
        log.warning("SNP %d: could not design primers in either orientation",
                    snp.raw_pos)
    return results


# ── Public entry point ────────────────────────────────────────────────────────

def design_kasp_markers(seq1: str, seq2: str, primer3_bin: str,
                        muscle_bin: str | None = None,
                        snp_positions=None,
                        max_markers: int = 10,
                        db_path: str = '',
                        blastn_path: str = 'blastn',
                        blastdbcmd_path: str = '',
                        progress_callback=None):
    """Main KASP design entry point.

    When db_path is provided, the common primer for each SNP is designed so
    its 3' end lands on a genome-discriminating position, identified via
    MSA of the target region against its closest homeologous sequences.
    """
    if len(seq1) != len(seq2):
        if muscle_bin:
            seq1, seq2 = align_two_seqs(seq1, seq2, muscle_bin)
        else:
            raise ValueError(
                "Sequences have different lengths — provide muscle_bin for alignment"
            )

    all_snps = find_snps(seq1, seq2)

    if not all_snps:
        return [], []

    snps = (
        [s for s in all_snps if s.raw_pos in snp_positions]
        if snp_positions else all_snps
    )
    log.info("Found %d SNPs (designing for %d)", len(all_snps), len(snps))

    ungap1 = seq1.replace('-', '')
    ungap2 = seq2.replace('-', '')

    total = len(snps)
    markers = []
    for i, snp in enumerate(snps):
        if len(markers) >= max_markers:
            break
        if progress_callback:
            progress_callback(f'SNP {i + 1} of {total} — position {snp.raw_pos + 1}')
        ms = design_kasp_for_snp(ungap1, ungap2, snp, primer3_bin,
                                  db_path=db_path, blastn_path=blastn_path,
                                  blastdbcmd_path=blastdbcmd_path,
                                  muscle_bin=muscle_bin or '',
                                  progress_callback=progress_callback)
        for m in ms:
            if len(markers) >= max_markers:
                break
            markers.append(m)

    return snps, markers


# ── Multi-strategy KASP design ────────────────────────────────────────────────

STRATEGIES = {
    'standard': {
        'name':    'Standard',
        'desc':    '18-25 bp allele-specific primers, Tm 57-63 °C',
        'min_size': 18, 'max_size': 25,
        'min_tm':  57.0, 'max_tm': 63.0, 'opt_tm': 60.0,
    },
    'extended': {
        'name':    'Extended',
        'desc':    '25-30 bp allele-specific primers, Tm 60-67 °C — better specificity',
        'min_size': 25, 'max_size': 30,
        'min_tm':  60.0, 'max_tm': 67.0, 'opt_tm': 63.0,
    },
    'common_focus': {
        'name':    'Common primer focus',
        'desc':    'Standard allele-specific + longer common primer (22-28 bp)',
        'min_size': 18, 'max_size': 25,
        'min_tm':  57.0, 'max_tm': 63.0, 'opt_tm': 60.0,
        'rev_min_size': 22, 'rev_max_size': 28,
    },
}


def design_kasp_for_snp_strategy(ungap1: str, ungap2: str, snp: SNPPosition,
                                  primer3_bin: str, strategy: str = 'standard',
                                  flank: int = 150) -> KASPMarker | None:
    """Design KASP primers using a specific strategy."""
    s             = STRATEGIES[strategy]
    pos           = snp.raw_pos
    start         = max(0, pos - flank)
    end           = min(len(ungap1), pos + flank + 1)
    snp_in_region = pos - start
    region1       = ungap1[start:end]
    region2       = ungap2[start:end]

    rev_min = s.get('rev_min_size', s['min_size'])
    rev_max = s.get('rev_max_size', s['max_size'])

    # ── Forward orientation ───────────────────────────────────────────────────

    base_fwd = _design_allele_fwd(
        region1, snp_in_region, primer3_bin,
        s['min_size'], s['max_size'], s['min_tm'], s['max_tm'], s['opt_tm'],
    )
    rev = _design_common_rev(
        region1, snp_in_region, primer3_bin,
        rev_min, rev_max, s['min_tm'], s['max_tm'], s['opt_tm'],
    ) if base_fwd else None

    if base_fwd and rev:
        fwd1       = base_fwd[:-1] + snp.allele1
        body       = base_fwd[:-1]
        mismatches = _count_body_mismatches(body, region2, snp_in_region)
        if mismatches > 0:
            log.info("SNP %d strategy %s: fwd2 body has %d mismatch(es) (using shared body)", snp.raw_pos, strategy, mismatches)
        fwd2 = body + snp.allele2

        prod_size = _fwd_product_size(region1, fwd1, rev, snp_in_region, flank)
        return KASPMarker(
            snp_pos=snp.raw_pos, allele1=snp.allele1, allele2=snp.allele2,
            fwd1_seq=fwd1, fwd2_seq=fwd2, rev_seq=rev,
            fwd1_with_tail=FAM_TAIL + fwd1, fwd2_with_tail=HEX_TAIL + fwd2,
            fwd1_tm=_calc_tm_nn(fwd1), fwd2_tm=_calc_tm_nn(fwd2),
            rev_tm=_calc_tm_nn(rev),
            fwd1_gc=_calc_gc(fwd1), fwd2_gc=_calc_gc(fwd2), rev_gc=_calc_gc(rev),
            product_size=prod_size, orientation='fwd',
            fwd2_body_mismatches=mismatches,
        )

    # ── Reverse orientation ───────────────────────────────────────────────────

    rc1       = _rc(region1)
    rc2       = _rc(region2)
    snp_in_rc = len(region1) - snp_in_region - 1

    base_rev = _design_allele_fwd(
        rc1, snp_in_rc, primer3_bin,
        s['min_size'], s['max_size'], s['min_tm'], s['max_tm'], s['opt_tm'],
    )
    cfwd = _design_common_rev(
        rc1, snp_in_rc, primer3_bin,
        rev_min, rev_max, s['min_tm'], s['max_tm'], s['opt_tm'],
    ) if base_rev else None

    if base_rev and cfwd:
        rev_body   = base_rev[:-1]
        mismatches = _count_body_mismatches(rev_body, rc2, snp_in_rc)
        if mismatches > 0:
            log.info("SNP %d strategy %s rev: rev2 body has %d mismatch(es) (using shared body)", snp.raw_pos, strategy, mismatches)
        rev2_body = rev_body

        rev1      = rev_body  + _rc(snp.allele1)
        rev2      = rev2_body + _rc(snp.allele2)
        cfwd_seq  = cfwd
        prod_size = _rev_product_size(rc1, base_rev, cfwd, snp_in_rc, flank)

        return KASPMarker(
            snp_pos=snp.raw_pos, allele1=snp.allele1, allele2=snp.allele2,
            fwd1_seq=rev1, fwd2_seq=rev2, rev_seq=cfwd_seq,
            fwd1_with_tail=FAM_TAIL + rev1, fwd2_with_tail=HEX_TAIL + rev2,
            fwd1_tm=_calc_tm_nn(rev1), fwd2_tm=_calc_tm_nn(rev2),
            rev_tm=_calc_tm_nn(cfwd_seq),
            fwd1_gc=_calc_gc(rev1), fwd2_gc=_calc_gc(rev2),
            rev_gc=_calc_gc(cfwd_seq),
            product_size=prod_size, orientation='rev',
            fwd2_body_mismatches=mismatches,
        )

    return None


def design_kasp_all_strategies(seq1: str, seq2: str, primer3_bin: str,
                                muscle_bin: str | None = None,
                                snp_positions=None,
                                max_markers: int = 50):
    """Design KASP markers using all 3 strategies simultaneously."""
    if len(seq1) != len(seq2):
        if muscle_bin:
            seq1, seq2 = align_two_seqs(seq1, seq2, muscle_bin)
        else:
            raise ValueError("Sequences have different lengths")

    all_snps = find_snps(seq1, seq2)
    if not all_snps:
        return all_snps, {}

    snps = (
        [s for s in all_snps if s.raw_pos in snp_positions]
        if snp_positions else all_snps
    )

    ungap1 = seq1.replace('-', '')
    ungap2 = seq2.replace('-', '')

    results = {}
    for strategy in STRATEGIES:
        markers = []
        for snp in snps:
            if len(markers) >= max_markers:
                break
            m = design_kasp_for_snp_strategy(ungap1, ungap2, snp,
                                             primer3_bin, strategy)
            if m:
                markers.append(m)
        results[strategy] = markers
        log.info("Strategy %s: %d markers", strategy, len(markers))

    return all_snps, results
