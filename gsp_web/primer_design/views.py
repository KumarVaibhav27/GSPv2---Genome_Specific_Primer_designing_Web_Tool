import os
import sys
import uuid
import time
import shutil
import threading as _threading
import subprocess
import tempfile
import json as _json
from pathlib import Path
import re as _re

# ── Async job state (file-based) ─────────────────────────────────────────────
# Each job's state is stored in <MEDIA_ROOT>/jobs/<job_id>/status.json so that
# any Gunicorn worker process can read or update it — no in-process dict needed.
# Atomic writes use a tmp file + os.replace() so readers never see partial JSON.

_JOB_TTL_SECONDS = 24 * 60 * 60  # 24 hours

def _job_dir(job_id: str) -> 'Path':
    return Path(settings.MEDIA_ROOT) / 'jobs' / job_id

def _read_job(job_id: str) -> 'dict | None':
    """Read job status from disk. Returns None if not found or unreadable."""
    try:
        return _json.loads((_job_dir(job_id) / 'status.json').read_text())
    except Exception:
        return None

def _write_job(job_id: str, data: dict) -> None:
    """Atomically write job status to disk (tmp + rename)."""
    d = _job_dir(job_id)
    d.mkdir(parents=True, exist_ok=True)
    tmp = d / 'status.tmp'
    tmp.write_text(_json.dumps(data))
    tmp.replace(d / 'status.json')

def _update_job(job_id: str, **kwargs) -> None:
    """Read-modify-write job status. Safe: one writer thread per job."""
    job = _read_job(job_id) or {}
    job.update(kwargs)
    _write_job(job_id, job)

def _cleanup_old_jobs():
    """Delete job folders older than TTL. Called lazily on each progress poll."""
    jobs_root = Path(settings.MEDIA_ROOT) / 'jobs'
    if not jobs_root.exists():
        return
    now = time.time()
    for job_dir in jobs_root.iterdir():
        if not job_dir.is_dir():
            continue
        job = _read_job(job_dir.name)
        if job and now - job.get('started_at', now) > _JOB_TTL_SECONDS:
            shutil.rmtree(job_dir, ignore_errors=True)

BLAST_STEPS = [
    (1, 'Reading & validating query'),
    (2, 'Running BLAST against genome'),
    (3, 'Filtering & checking hits'),
    (4, 'Extracting homeolog sequences'),
    (5, 'Building multiple sequence alignment'),
    (6, 'Designing primer candidates (Primer3)'),
    (7, 'Genome-wide BLAST specificity filter'),
    (8, 'Finalising results'),
]

KASP_STEPS = [
    (1, 'Aligning sequences'),
    (2, 'Detecting SNPs'),
    (3, 'Designing primers'),
    (4, 'Genome BLAST specificity check'),
    (5, 'Off-target amplicon & thermodynamics'),
    (6, 'Finalising results'),
]

def _progress(job_id: str, step: int, detail: str = '') -> None:
    """Update progress for a running job.  Called from background thread."""
    step_name = next((n for s, n in BLAST_STEPS if s == step), f'Step {step}')
    _update_job(job_id, step=step, step_name=step_name, detail=detail)

def _kasp_progress(job_id: str, step: int, detail: str = '') -> None:
    """Update progress for a running KASP job.  Called from background thread."""
    step_name = next((n for s, n in KASP_STEPS if s == step), f'Step {step}')
    _update_job(job_id, step=step, step_name=step_name, detail=detail)

def _kasp_detail(job_id: str, detail: str) -> None:
    """Update only the detail text for a KASP job without advancing the step number."""
    _update_job(job_id, detail=detail)

def calc_thermo(seq1, seq2=None, mode='END1', ntthal_bin=None):
    if ntthal_bin is None:
        ntthal_bin = getattr(settings, 'NTTHAL_PATH', '/opt/homebrew/bin/ntthal')
    import subprocess
    cmd = [ntthal_bin, '-a', mode, '-s1', seq1]
    if seq2:
        cmd += ['-s2', seq2]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        out = r.stdout + r.stderr
        dg = _re.search(r'dG\s*=\s*([-\d.]+)', out)
        tm = _re.search(r't\s*=\s*([-\d.]+)', out)
        dg_val = round(float(dg.group(1)) / 1000, 2) if dg else None
        tm_val = round(float(tm.group(1)), 1) if tm else None
        return {'dg': dg_val, 'tm': tm_val}
    except:
        return {'dg': None, 'tm': None}

def _thermo_risk_label(dg):
    """Return dG value with risk label e.g. '-7.2 (Moderate)'. Returns 'None detected' if dg is None."""
    if dg is None:
        return 'None detected'
    label = 'Low risk' if dg > -6 else ('Moderate' if dg > -9 else 'High risk ⚠')
    return f'{dg} ({label})'

def _thermo_summary(dg_values):
    """Return 'X/Y' checks-passed string from a list of dG values (None allowed)."""
    non_none = [v for v in dg_values if v is not None]
    passed = sum(1 for v in non_none if v > -6)
    return f'{passed}/{len(non_none)}'

def _autofit_xlsx(wb, padding=2, max_width=60):
    """Set each column width to its longest cell value to avoid overlap."""
    for ws in wb.worksheets:
        for col in ws.columns:
            max_len = max(
                (len(str(cell.value)) for cell in col if cell.value is not None),
                default=8
            )
            ws.column_dimensions[col[0].column_letter].width = min(max_len + padding, max_width)

from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse, JsonResponse
from django.conf import settings

# Add gsp_python to path
sys.path.insert(0, str(Path(settings.GSP_PYTHON_DIR)))


def _validate_fasta(text: str, min_length: int = 100) -> str | None:
    """
    Validate a FASTA string before it is passed to MUSCLE or Primer3.
    Returns None if valid, or a human-readable error string if not.

    Checks:
    - At least one sequence is present
    - Each sequence is at least min_length bp
    - Sequences contain only valid IUPAC DNA characters (ACGTNRYMKSWHBVD and gaps)
    - No RNA (U bases) — common copy-paste mistake
    """
    valid_bases = set('ACGTNRYMKSWHBVDacgtnrymkswhbvd-')
    lines   = text.strip().splitlines()
    headers = 0
    current_seq_chars: list[str] = []
    seq_lengths: list[int] = []

    def _flush():
        if current_seq_chars:
            seq_lengths.append(len(current_seq_chars))
        current_seq_chars.clear()

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            _flush()
            headers += 1
            current_seq_chars = []
        else:
            for ch in line:
                if ch.upper() == 'U':
                    return ('Sequence contains RNA base "U". '
                            'Please convert RNA to DNA (replace U with T) before submitting.')
                if ch not in valid_bases:
                    return (f'Invalid character "{ch}" found in sequence. '
                            f'Only standard IUPAC DNA bases (A, C, G, T, N, and ambiguity codes) are allowed.')
                current_seq_chars.append(ch)
    _flush()

    if headers == 0:
        # Treat the whole input as a bare sequence (no header line)
        seq_lengths = [sum(1 for c in text if c in valid_bases)]
        if not seq_lengths[0]:
            return 'No sequence found. Please provide a DNA sequence.'

    if not seq_lengths:
        return 'No sequence found. Please provide a DNA sequence.'

    for i, length in enumerate(seq_lengths, 1):
        if length < min_length:
            return (f'Sequence {i} is too short ({length} bp). '
                    f'A minimum of {min_length} bp is required for primer design.')

    return None   # valid


def _normalize_chrom(name: str) -> str:
    """
    Return the shortest unique identifier for a chromosome name.
    Strips NCBI-style pipe-delimited prefixes and colon suffixes so that
    '5A', 'chr5A', 'lcl|5A', 'gb|5A|' and '5A:1-700000000' all normalise
    to '5A'.
    """
    # Strip everything after the first colon (range suffixes)
    name = name.split(':')[0]
    # Strip pipe-delimited accession prefixes (NCBI format)
    parts = name.split('|')
    name = parts[-1] if parts[-1] else parts[-2]
    return name.strip()


def _build_target_chroms(query_chrom: str) -> set:
    """
    Build a set of all name variants for the query chromosome so that the
    genome-wide BLAST filter can match hits regardless of naming convention.
    Replaces the old ad-hoc 4-variant construction in run_pipeline.
    """
    variants = set()
    variants.add(query_chrom)
    variants.add(_normalize_chrom(query_chrom))
    # Also add without version suffixes like .1, .2
    base = _normalize_chrom(query_chrom)
    if '.' in base:
        variants.add(base.rsplit('.', 1)[0])
    return variants


def index(request):
    """Home page — sequence input form."""
    databases = list(settings.BLAST_DATABASES.keys())
    return render(request, 'primer_design/index.html', {
        'databases': databases,
        'databases_grouped': settings.BLAST_DATABASES_GROUPED,
    })



def _primer_3end_locations(primer_seq, db_path, end_len=12, max_locs=10):
    """Check if the last end_len bases of primer hit multiple genome locations.
    Returns list of locations where 3' end matches at 100% identity."""
    if not db_path or len(primer_seq) < end_len:
        return []
    end_seq = primer_seq[-end_len:]
    try:
        import subprocess as _sp
        r = _sp.run([
            settings.BLASTN_PATH, '-db', db_path,
            '-query', '-', '-task', 'blastn-short',
            '-word_size', '4', '-outfmt', '6',
            '-max_target_seqs', '50', '-evalue', '1000',
            '-dust', 'no',
        ], input=f'>p\n{end_seq}\n',
        capture_output=True, text=True, timeout=30)
        locs, seen = [], set()
        for hit in r.stdout.strip().split('\n'):
            if not hit: continue
            parts = hit.split('\t')
            if len(parts) >= 9:
                pident  = float(parts[2])
                aln_len = int(parts[3])
                if pident >= 90.0 and aln_len >=8:
                    chrom = parts[1].split(':')[0]
                    if chrom not in seen:
                        seen.add(chrom)
                        locs.append({'chrom': chrom, 'start': int(parts[8]),
                                     'end': int(parts[9]), 'identity': pident})
                    if len(locs) >= max_locs: break
        return locs
    except:
        return []

def _check_pair_off_target(left_seq, right_seq, db_path, max_amplicon=1000):
    """Check if primer pair can amplify off-target regions < max_amplicon bp.
    Returns True if pair is safe (no off-target amplification).
    Returns False on exception — conservative: if we cannot check, assume unsafe.
    """
    if not db_path:
        return True
    try:
        left_locs  = _primer_locations(left_seq,  db_path, max_locs=20)
        right_locs = _primer_locations(right_seq, db_path, max_locs=20)
        # Check if any left+right pair is within max_amplicon on same chromosome
        for ll in left_locs:
            for rl in right_locs:
                if ll['chrom'] == rl['chrom']:
                    dist = abs(int(rl['start']) - int(ll['start']))
                    if 0 < dist < max_amplicon:
                        return False  # off-target amplification possible
        return True
    except Exception:
        # If the BLAST check itself fails we cannot confirm safety — return False
        return False

def _check_primer_specificity(primer_seq, db_path):
    """
    Full specificity check per specification.
    Returns (is_3end_specific, is_full_specific, locations).

    - is_full_specific: full primer hits EXACTLY 1 genome location (the target).
      >1 hit means the primer is not genome-specific.
    - is_3end_specific: the last 12bp hit EXACTLY 1 genome location at 100%.
      This is the stricter, more biologically meaningful criterion because
      the 3' end determines where polymerase initiates extension.
    """
    if not db_path:
        return True, True, []
    full_locs = _primer_locations(primer_seq, db_path, max_locs=10)
    end_locs  = _primer_3end_locations(primer_seq, db_path)
    # A primer is specific only if it maps to EXACTLY 1 location genome-wide.
    is_full  = len(full_locs) == 1
    is_3end  = len(end_locs)  == 1
    return is_3end, is_full, full_locs

def _primer_locations(primer_seq, db_path, max_locs=5):
    if not db_path:
        return []
    try:
        import subprocess as _sp
        r = _sp.run([
            settings.BLASTN_PATH, '-db', db_path,
            '-query', '-', '-task', 'blastn-short',
            '-word_size', '11', '-outfmt', '6',
            '-max_target_seqs', '100', '-evalue', '1000',
            '-dust', 'no',
        ], input=f'>p\n{primer_seq}\n',
        capture_output=True, text=True, timeout=30)
        locs, seen = [], set()
        for hit in r.stdout.strip().split('\n'):
            if not hit: continue
            parts = hit.split('\t')
            if len(parts) >= 9 and float(parts[2]) >= 95 and int(parts[3]) >= int(0.85 * len(primer_seq)):
                chrom = parts[1].split(':')[0]
                if chrom not in seen:
                    seen.add(chrom)
                    locs.append({'chrom': chrom, 'start': int(parts[8]),
                                 'end': int(parts[9]), 'identity': float(parts[2])})
                if len(locs) >= max_locs: break
        return locs
    except:
        return []

def _compute_uniqueness_map(seq: str, db_path: str) -> list:
    """
    BLAST seq against the genome. Return per-position boolean: True = unique
    (maps to ≤1 chromosome), False = repetitive/multi-copy.
    """
    if not db_path or not seq:
        return [True] * len(seq)
    try:
        import subprocess as _sp
        r = _sp.run([
            settings.BLASTN_PATH, '-db', db_path,
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


# ── Background BLAST pipeline ─────────────────────────────────────────────────

def _set_job_error(job_id: str, msg: str) -> None:
    _update_job(job_id, status='error', error=msg)

def _set_job_done(job_id: str, result: dict) -> None:
    _update_job(job_id, status='done', result=result)


def _run_blast_pipeline(job_id: str, workdir: Path,
                        query_path: Path, db_path: str, db_name: str,
                        amp_min: int, amp_max: int) -> None:
    """
    Long-running BLAST + primer design pipeline.
    Runs in a background thread.  Writes progress and final result to
    status.json in the job workdir for the polling endpoint.
    """
    import sys as _sys
    _sys.path.insert(0, str(Path(settings.GSP_PYTHON_DIR)))

    try:
        # ── Step 1: read query ────────────────────────────────────────────────
        _progress(job_id, 1, 'Reading query sequence')
        query_seq  = query_path.read_text()
        query_len  = sum(len(l.strip()) for l in query_seq.split('\n')
                         if l.strip() and not l.startswith('>'))
        if query_len < 100:
            _set_job_error(job_id,
                f'Query sequence too short ({query_len} bp). Minimum 100 bp required.')
            return

        # ── Step 2: run BLAST ─────────────────────────────────────────────────
        _progress(job_id, 2, f'BLASTing {query_len:,} bp against {db_name}')
        blast_out = workdir / 'blast.tab'
        blast_cmd = [
            settings.BLASTN_PATH,
            '-query', str(query_path), '-db', db_path,
            '-out',   str(blast_out),
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
            '-evalue', '1e-5', '-dust', 'yes', '-soft_masking', 'true',
            '-perc_identity', '80', '-num_alignments', '20', '-num_threads', '4',
        ]
        r = subprocess.run(blast_cmd, capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            _set_job_error(job_id, f'BLAST failed: {r.stderr[:300]}')
            return
        if not blast_out.exists() or blast_out.stat().st_size == 0:
            _set_job_error(job_id, 'BLAST returned no hits. Try a different sequence or database.')
            return

        # ── Step 3: filter & check ────────────────────────────────────────────
        _progress(job_id, 3, 'Filtering hits and checking for repeats')
        min_aln_len = min(max(200, int(query_len * 0.05)), 600)
        filtered_blast = workdir / 'blast_filtered.tab'
        with open(blast_out) as fin, open(filtered_blast, 'w') as fout:
            for line in fin:
                parts = line.strip().split('\t')
                if len(parts) >= 12 and int(parts[3]) >= min_aln_len:
                    fout.write(line)
        if not filtered_blast.exists() or filtered_blast.stat().st_size == 0:
            _set_job_error(job_id,
                f'No significant BLAST hits (required alignment ≥ {min_aln_len} bp). '
                f'Try a longer query sequence.')
            return

        # Repetitiveness check
        from collections import Counter as _Ctr
        _full_hit_floor = max(300, int(query_len * 0.5))
        _full_hits = _Ctr()
        with open(filtered_blast) as _fh:
            for _ln in _fh:
                _p = _ln.strip().split('\t')
                if len(_p) >= 4 and int(_p[3]) >= _full_hit_floor:
                    _full_hits[_p[1]] += 1
        _repeat_warning = None
        if _full_hits:
            _top_chrom, _top_count = _full_hits.most_common(1)[0]
            if _top_count > 5:
                _short = _top_chrom.split('|')[-1]
                _repeat_warning = (
                    f'Warning: this region appears repetitive ({_short} has '
                    f'{_top_count} near-full-length hits). Genome-specific '
                    f'primers may not be achievable — proceeding anyway.'
                )
                _progress(job_id, 3, _repeat_warning)

        # ── Step 4: extract homeolog sequences ────────────────────────────────
        from blast     import parse_blast_tab
        from alignment import _calc_bed_coords
        groups = parse_blast_tab(str(filtered_blast))
        if not groups:
            _set_job_error(job_id, 'No valid BLAST groups found.')
            return

        # Use all detected chromosomes (up to 30 to cap runtime for extreme cases).
        # This ensures homeologs AND any other similar loci are all included in
        # the MSA — more discriminating positions found = better primer specificity.
        MAX_HOMEOLOGS = 30
        _group_raw  = groups[0]
        _n_detected = _group_raw.num_chroms
        _n_use      = min(_n_detected, MAX_HOMEOLOGS)
        group       = _group_raw.select_top_n(_n_use)

        _flank = 15
        _progress(job_id, 4,
                  f'Extracting {_n_use} homeolog sequence{"s" if _n_use > 1 else ""} '
                  f'({_n_detected} chromosomes hit)')

        extracted   = {}
        blastdbcmd  = settings.BLASTDBCMD_PATH
        genome_fasta = settings.BLAST_GENOME_FASTAS.get(db_name)
        use_fasta_fallback = genome_fasta and Path(genome_fasta).exists()
        if use_fasta_fallback:
            from alignment import read_fasta as _rfa
            _all_seqs_raw   = _rfa(genome_fasta)
            _short_name_map = {h.split()[0]: s for h, s in _all_seqs_raw.items()}

        for chrom, hits in group.hits_by_chrom.items():
            # Use only the best (longest alignment) hit for coordinate
            # calculation.  Secondary partial hits on the same chromosome can
            # be thousands of bp away, inflating the extraction window to
            # include intronic / unrelated sequence and producing primers that
            # fall outside the user's pasted query.
            best_hit = max(hits, key=lambda h: abs(h.send - h.sstart))
            start, end = _calc_bed_coords([best_hit], _flank, 10**10)
            strand = hits[0].strand
            if use_fasta_fallback:
                if chrom not in _short_name_map:
                    continue
                sub = _short_name_map[chrom][start-1:end]
                if strand == '-':
                    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
                    sub  = sub.translate(comp)[::-1]
                extracted[chrom] = sub
            else:
                strand_flag = 'minus' if strand == '-' else 'plus'
                cmd = [blastdbcmd, '-db', db_path, '-entry', chrom,
                       '-range', f'{start}-{end}', '-strand', strand_flag, '-outfmt', '%f']
                rc = subprocess.run(cmd, capture_output=True, text=True)
                if rc.returncode != 0 or not rc.stdout.strip():
                    continue  # skip; other chroms may still work
                seq = ''.join(l for l in rc.stdout.strip().split('\n') if not l.startswith('>'))
                extracted[chrom] = seq

        # ── Single-copy mode ──────────────────────────────────────────────────
        if len(extracted) == 1:
            from primers import design_primers_single_copy
            _sc_chrom  = list(extracted.keys())[0]
            _sc_seq    = list(extracted.values())[0]
            _tc = {_normalize_chrom(n) for n in extracted.keys()}
            _tc.update(extracted.keys())
            from primers import DEFAULT_SETTINGS as _P3_DEFAULTS
            _sc_settings = dict(_P3_DEFAULTS)
            _sc_settings.update({
                'PRIMER_PRODUCT_SIZE_RANGE': f'{amp_min}-{amp_max}',
                'PRIMER_MIN_GC': '35.0', 'PRIMER_MAX_GC': '65.0',
                'PRIMER_MIN_TM': '50.0', 'PRIMER_OPT_TM': '60.0',
                'PRIMER_MAX_TM': '65.0', 'PRIMER_PAIR_MAX_DIFF_TM': '2.0',
                'PRIMER_MAX_POLY_X': '3',
                'PRIMER_MIN_SIZE': '18',  'PRIMER_MAX_SIZE': '25',
            })
            _progress(job_id, 4, 'Single-copy gene — MSA alignment skipped')
            _progress(job_id, 5, '')
            _progress(job_id, 6, f'Designing primers for {_sc_chrom[:40]}')
            _progress(job_id, 7, 'Genome-wide BLAST specificity filter (may take several minutes)')
            pairs = design_primers_single_copy(
                sequence=_sc_seq, seq_name=_sc_chrom,
                primer3_bin=settings.PRIMER3_PATH,
                settings=_sc_settings, blast_db=db_path,
                blastn_path=settings.BLASTN_PATH, target_chroms=_tc,
                require_3end_unique=True,
                included_region=(15, max(1, len(_sc_seq) - 30)),
            )
            # ── Query-sequence filter (single-copy mode) ──────────────────────
            _sc_query_dna = ''.join(
                l.strip().upper() for l in query_seq.split('\n')
                if l.strip() and not l.startswith('>')
            )
            def _sc_primer_in_query(pair):
                from primers import reverse_complement as _rc
                return (pair.left_seq.upper() in _sc_query_dna and
                        _rc(pair.right_seq).upper() in _sc_query_dna)
            _sc_before = len(pairs)
            pairs = [p for p in pairs if _sc_primer_in_query(p)]
            if _sc_before - len(pairs):
                import logging as _log
                _log.getLogger(__name__).info(
                    'Single-copy query-sequence filter removed %d pair(s) '
                    'outside the pasted query', _sc_before - len(pairs)
                )

            if not pairs:
                _set_job_done(job_id, {
                    'success': True, 'mode': 'single_copy', 'db_name': db_name,
                    'total_pairs': 0,
                    'subgenomes': [{'name': _sc_chrom, 'count': 0, 'pairs': []}],
                    'top_pairs': [],
                    'msa_data': None,
                    'job_id': job_id,
                    'warning': (
                        'No genome-specific primer pairs found for this region. '
                        'The sequence may be too conserved or repetitive for '
                        'genome-specific primer design. '
                        'Try adjusting the amplicon size range.'
                    ),
                })
                return
            # Cap to top 20 pairs before writing any output
            pairs = pairs[:20]

            _progress(job_id, 8, f'Writing output ({len(pairs)} pairs)')
            import openpyxl as _ox, csv as _csv_mod
            xlsx_path = workdir / 'primers.xlsx'
            _wb = _ox.Workbook(); _wb.remove(_wb.active)
            _hdrs = ['Rank','Sequence','Left primer (5->3)','Left start','Left len',
                     'Left Tm','Left GC%','Right primer (5->3)','Right start','Right len',
                     'Right Tm','Right GC%','Product size','Off-target amplicon',
                     'Heterodimer dG','L Hairpin dG','R Hairpin dG',
                     'L Homodimer dG','R Homodimer dG','Thermo checks passed','Score']
            _ws = _wb.create_sheet(title=_sc_chrom[:30]); _ws.append(_hdrs)
            _csv_rows_sc = []
            for _rk, _p in enumerate(pairs, 1):
                _th_sc = {
                    'heterodimer':    calc_thermo(_p.left_seq, _p.right_seq, 'END1').get('dg'),
                    'left_hairpin':   calc_thermo(_p.left_seq, mode='HAIRPIN').get('dg'),
                    'right_hairpin':  calc_thermo(_p.right_seq, mode='HAIRPIN').get('dg'),
                    'left_homodimer': calc_thermo(_p.left_seq, _p.left_seq, 'ANY').get('dg'),
                    'right_homodimer':calc_thermo(_p.right_seq, _p.right_seq, 'ANY').get('dg'),
                }
                _ot_sc = 'Yes' if getattr(_p, 'off_target_amplicon', False) else 'No'
                _ws.append([_rk, _p.seq_name, _p.left_seq, _p.left_start, _p.left_len,
                             float(_p.left_tm), float(_p.left_gc),
                             _p.right_seq, _p.right_start, _p.right_len,
                             float(_p.right_tm), float(_p.right_gc),
                             _p.product_size, _ot_sc,
                             _thermo_risk_label(_th_sc['heterodimer']),
                             _thermo_risk_label(_th_sc['left_hairpin']),
                             _thermo_risk_label(_th_sc['right_hairpin']),
                             _thermo_risk_label(_th_sc['left_homodimer']),
                             _thermo_risk_label(_th_sc['right_homodimer']),
                             _thermo_summary([_th_sc['heterodimer'], _th_sc['left_hairpin'],
                                              _th_sc['right_hairpin'], _th_sc['left_homodimer'],
                                              _th_sc['right_homodimer']]),
                             round(_p.score, 2)])
                _csv_rows_sc.append([_rk, _p.seq_name, _p.left_seq, _p.left_start, _p.left_len,
                                     float(_p.left_tm), float(_p.left_gc),
                                     _p.right_seq, _p.right_start, _p.right_len,
                                     float(_p.right_tm), float(_p.right_gc), _p.product_size, _ot_sc,
                                     _thermo_risk_label(_th_sc['heterodimer']),
                                     _thermo_risk_label(_th_sc['left_hairpin']),
                                     _thermo_risk_label(_th_sc['right_hairpin']),
                                     _thermo_risk_label(_th_sc['left_homodimer']),
                                     _thermo_risk_label(_th_sc['right_homodimer']),
                                     _thermo_summary([_th_sc['heterodimer'], _th_sc['left_hairpin'],
                                                      _th_sc['right_hairpin'], _th_sc['left_homodimer'],
                                                      _th_sc['right_homodimer']]),
                                     round(_p.score, 2)])
            _autofit_xlsx(_wb)
            _wb.save(str(xlsx_path))

            # Write CSV — same format as xlsx but without 3end columns
            _csv_path = workdir / 'primers.csv'
            _csv_hdrs_sc = ['Rank','Sequence','Left primer (5->3)','Left start','Left len',
                            'Left Tm','Left GC%','Right primer (5->3)','Right start','Right len',
                            'Right Tm','Right GC%','Product size','Off-target amplicon',
                            'Heterodimer dG','L Hairpin dG','R Hairpin dG',
                            'L Homodimer dG','R Homodimer dG','Thermo checks passed','Score']
            with open(str(_csv_path), 'w', newline='') as _cf:
                _cw = _csv_mod.writer(_cf)
                _cw.writerow(_csv_hdrs_sc)
                _cw.writerows(_csv_rows_sc)

            # Build subgenomes + top_pairs in the SAME format as MSA mode so the
            # frontend can render the primer table without any special-casing.
            _sc_pairs_out = []
            for _rk2, _p in enumerate(pairs, 1):
                _sc_pairs_out.append({
                    'rank': _rk2, 'seq_name': _p.seq_name,
                    'left_seq': _p.left_seq, 'left_start': _p.left_start,
                    'left_len': _p.left_len, 'left_tm': float(_p.left_tm),
                    'left_gc': float(_p.left_gc),
                    'right_seq': _p.right_seq, 'right_start': _p.right_start,
                    'right_len': _p.right_len, 'right_tm': float(_p.right_tm),
                    'right_gc': float(_p.right_gc),
                    'product_size': _p.product_size,
                    # Frontend uses 'left_diff'/'right_diff', not 'left_diff_count'
                    'left_diff':  _p.left_diff_count,
                    'right_diff': _p.right_diff_count,
                    'left_unique_pos':  _p.left_unique_pos  or [],
                    'right_unique_pos': _p.right_unique_pos or [],
                    # Frontend uses 'left_3end'/'right_3end'
                    'left_3end':  _p.left_end_genome_specific,
                    'right_3end': _p.right_end_genome_specific,
                    'left_genome_specific':  _p.left_genome_specific,
                    'right_genome_specific': _p.right_genome_specific,
                    'left_genome_hits':  _p.left_genome_hits,
                    'right_genome_hits': _p.right_genome_hits,
                    'off_target_amplicon': getattr(_p, 'off_target_amplicon', False),
                    'score': round(_p.score, 2),
                    'thermo': {
                        'heterodimer':    calc_thermo(_p.left_seq, _p.right_seq, 'END1'),
                        'left_hairpin':   calc_thermo(_p.left_seq, mode='HAIRPIN'),
                        'right_hairpin':  calc_thermo(_p.right_seq, mode='HAIRPIN'),
                        'left_homodimer': calc_thermo(_p.left_seq, _p.left_seq, 'ANY'),
                        'right_homodimer':calc_thermo(_p.right_seq, _p.right_seq, 'ANY'),
                    },
                    'left_locations': [], 'right_locations': [],
                })
            _set_job_done(job_id, {
                'success': True, 'mode': 'single_copy', 'db_name': db_name,
                'total_pairs': len(pairs),
                'subgenomes': [{'name': _sc_chrom, 'count': len(pairs),
                                'pairs': _sc_pairs_out}],
                'top_pairs': [{'rank': _p['rank'], 'seq_name': _p['seq_name'],
                               'left_seq': _p['left_seq'], 'left_start': _p['left_start'],
                               'left_len': _p['left_len'], 'left_tm': _p['left_tm'],
                               'left_gc': _p['left_gc'],
                               'right_seq': _p['right_seq'], 'right_start': _p['right_start'],
                               'right_len': _p['right_len'], 'right_tm': _p['right_tm'],
                               'right_gc': _p['right_gc'],
                               'product_size': _p['product_size'],
                               'left_diff': _p['left_diff'], 'right_diff': _p['right_diff'],
                               'score': _p['score']}
                              for _p in _sc_pairs_out[:50]],
                'msa_data': None,   # no alignment in single-copy mode
                'job_id': job_id,
            })
            return

        if len(extracted) == 0:
            _set_job_error(job_id,
                'No sequences could be extracted. '
                'Check that chromosome names in BLAST output match the database.')
            return

        fasta_path  = workdir / 'sequences.fas'
        query_chrom = list(extracted.keys())[0]
        with open(fasta_path, 'w') as _f:
            for _name, _seq in extracted.items():
                _f.write(f'>{_name}\n{_seq}\n')

        # ── Step 5: multiple sequence alignment ───────────────────────────────
        from alignment import build_msa
        _progress(job_id, 5,
                  f'Aligning {len(extracted)} sequences with MUSCLE '
                  f'({"up to " if len(extracted)>2 else ""}{sum(len(s) for s in extracted.values())//len(extracted):,} bp avg)')
        msa, unique = build_msa(None, str(fasta_path), settings.MUSCLE_PATH)
        if not unique:
            _set_job_error(job_id,
                'No unique positions found between homeologs — sequences are too similar.')
            return

        # ── Step 6+7: Primer3 + genome BLAST filter ───────────────────────────
        amp_min = max(50,  min(amp_min, 5000))
        amp_max = max(amp_min + 50, min(amp_max, 5000))
        from primers import design_primers, DEFAULT_SETTINGS
        primer_settings = dict(DEFAULT_SETTINGS)
        primer_settings['PRIMER_PRODUCT_SIZE_RANGE'] = f'{amp_min}-{amp_max}'
        primer_settings['PRIMER_MIN_GC']           = '35.0'
        primer_settings['PRIMER_MAX_GC']           = '65.0'
        primer_settings['PRIMER_MIN_TM']           = '50.0'
        primer_settings['PRIMER_OPT_TM']           = '60.0'
        primer_settings['PRIMER_MAX_TM']           = '65.0'
        primer_settings['PRIMER_PAIR_MAX_DIFF_TM'] = '2.0'
        primer_settings['PRIMER_MAX_POLY_X']       = '3'
        primer_settings['PRIMER_MIN_SIZE']         = '18'
        primer_settings['PRIMER_MAX_SIZE']         = '25'

        from alignment import read_fasta as _rf2
        _full_seqs = _rf2(str(fasta_path))

        # Build target_chroms: include all extracted chromosomes with both raw and
        # normalised names so the genome-wide BLAST filter doesn't reject homeolog hits.
        # Normalise as split(':')[0] first (strip range suffixes), then split('|')[-1]
        # (strip NCBI prefix) — same order as _normalize_chrom().
        _target_chroms = set()
        for _n in extracted.keys():
            _target_chroms.add(_n)
            _target_chroms.add(_normalize_chrom(_n))

        # NOTE: design_primers runs Primer3 (step 6) then genome BLAST filter (step 7)
        # internally with no callback hooks.  Pre-announce step 7 so the UI shows the
        # correct label for the long genome-wide BLAST wait.
        _progress(job_id, 6, f'Primer3: {len(unique)} discriminating positions across {msa.num_seqs} sequences')
        _progress(job_id, 7,
                  f'Genome-wide BLAST filter against {db_name} '
                  f'(this step typically takes several minutes)')

        pairs = design_primers(
            msa, unique, settings.PRIMER3_PATH,
            settings            = primer_settings,
            full_seqs           = _full_seqs,
            blast_db            = db_path,
            blastn_path         = settings.BLASTN_PATH,
            target_chroms       = _target_chroms,
            require_3end_unique = True,
            primary_target      = query_chrom,
        )
        if not pairs:
            _set_job_done(job_id, {
                'success': True, 'mode': 'blast', 'db_name': db_name,
                'total_pairs': 0, 'subgenomes': [], 'top_pairs': [],
                'msa_data': None, 'job_id': job_id,
                'warning': (
                    'No genome-specific primer pairs found for this region. '
                    'The sequence may be too conserved or repetitive for '
                    'genome-specific primer design.'
                ),
            })
            return

        pairs = [p for p in pairs if p.seq_name == query_chrom]

        # ── Query-sequence filter ─────────────────────────────────────────────
        # Primers must be substrings of the user's pasted query.
        # This catches any primers that Primer3 placed outside the query region
        # (e.g. due to secondary BLAST hits inflating the extraction window).
        _query_dna = ''.join(
            l.strip().upper() for l in query_seq.split('\n')
            if l.strip() and not l.startswith('>')
        )
        def _primer_in_query(pair):
            from primers import reverse_complement as _rc
            l_fwd = pair.left_seq.upper()  in _query_dna
            r_rc  = _rc(pair.right_seq).upper() in _query_dna
            return l_fwd and r_rc
        _before = len(pairs)
        pairs = [p for p in pairs if _primer_in_query(p)]
        if len(pairs) < _before:
            import logging as _log
            _log.getLogger(__name__).info(
                'Query-sequence filter removed %d pair(s) whose primers were '
                'outside the pasted query (likely from secondary BLAST hits)',
                _before - len(pairs)
            )

        if not pairs:
            _set_job_done(job_id, {
                'success': True, 'mode': 'blast', 'db_name': db_name,
                'total_pairs': 0, 'subgenomes': [], 'top_pairs': [],
                'msa_data': None, 'job_id': job_id,
                'warning': (
                    f'No primers found specifically for the query chromosome '
                    f'({query_chrom}). The region may be too similar to other '
                    f'homeologs for genome-specific primer design. '
                    f'Try adjusting the amplicon size range.'
                ),
            })
            return

        # Cap to top 20 pairs
        pairs = pairs[:20]

        # ── Step 8: build output ──────────────────────────────────────────────
        _progress(job_id, 8, f'Writing Excel output ({len(pairs)} pairs)')
        import openpyxl
        from collections import defaultdict as _dd
        from output import write_csv
        import re as _re3

        xlsx_path = workdir / 'primers.xlsx'
        wb = openpyxl.Workbook(); wb.remove(wb.active)
        headers = ['Rank','Sequence','Left primer (5->3)','Left start','Left len',
                   'Left Tm','Left GC%','Right primer (5->3)','Right start','Right len',
                   'Right Tm','Right GC%','Product size','Off-target amplicon',
                   'Heterodimer dG','L Hairpin dG','R Hairpin dG',
                   'L Homodimer dG','R Homodimer dG','Thermo checks passed','Score']
        _by_seq = _dd(list)
        for _p in pairs:
            _by_seq[_p.seq_name].append(_p)
        import csv as _csv_mod2
        _csv_rows_mc = []
        for _sn, _sps in _by_seq.items():
            _cm = _re3.search(r'([1-7][ABDSabds])', _sn)
            _sht = _cm.group(1).upper() if _cm else _re3.sub(r'[:/\\?*\[\]]', '_', _sn.split()[0])[:31]
            _ws2 = wb.create_sheet(title=_sht); _ws2.append(headers)
            for _p in _sps:
                _th = {
                    'heterodimer':    calc_thermo(_p.left_seq, _p.right_seq, 'END1').get('dg'),
                    'left_hairpin':   calc_thermo(_p.left_seq, mode='HAIRPIN').get('dg'),
                    'right_hairpin':  calc_thermo(_p.right_seq, mode='HAIRPIN').get('dg'),
                    'left_homodimer': calc_thermo(_p.left_seq, _p.left_seq, 'ANY').get('dg'),
                    'right_homodimer':calc_thermo(_p.right_seq, _p.right_seq, 'ANY').get('dg'),
                }
                _rk_mc = _p.rank if hasattr(_p, 'rank') else _sps.index(_p)+1
                _ot_mc = 'Yes' if getattr(_p, 'off_target_amplicon', False) else 'No'
                _ws2.append([_rk_mc,
                              _p.seq_name, _p.left_seq, _p.left_start, _p.left_len,
                              _p.left_tm, _p.left_gc, _p.right_seq, _p.right_start,
                              _p.right_len, _p.right_tm, _p.right_gc, _p.product_size, _ot_mc,
                              _thermo_risk_label(_th['heterodimer']),
                              _thermo_risk_label(_th['left_hairpin']),
                              _thermo_risk_label(_th['right_hairpin']),
                              _thermo_risk_label(_th['left_homodimer']),
                              _thermo_risk_label(_th['right_homodimer']),
                              _thermo_summary([_th['heterodimer'], _th['left_hairpin'],
                                               _th['right_hairpin'], _th['left_homodimer'],
                                               _th['right_homodimer']]),
                              round(_p.score, 2)])
                _csv_rows_mc.append([_rk_mc,
                                     _p.seq_name, _p.left_seq, _p.left_start, _p.left_len,
                                     _p.left_tm, _p.left_gc, _p.right_seq, _p.right_start,
                                     _p.right_len, _p.right_tm, _p.right_gc, _p.product_size, _ot_mc,
                                     _thermo_risk_label(_th['heterodimer']),
                                     _thermo_risk_label(_th['left_hairpin']),
                                     _thermo_risk_label(_th['right_hairpin']),
                                     _thermo_risk_label(_th['left_homodimer']),
                                     _thermo_risk_label(_th['right_homodimer']),
                                     _thermo_summary([_th['heterodimer'], _th['left_hairpin'],
                                                      _th['right_hairpin'], _th['left_homodimer'],
                                                      _th['right_homodimer']]),
                                     round(_p.score, 2)])
        _autofit_xlsx(wb)
        wb.save(str(xlsx_path))

        # Write CSV — same format as xlsx but without 3end columns
        csv_path = workdir / 'primers.csv'
        _csv_hdrs_mc = ['Rank','Sequence','Left primer (5->3)','Left start','Left len',
                        'Left Tm','Left GC%','Right primer (5->3)','Right start','Right len',
                        'Right Tm','Right GC%','Product size','Off-target amplicon',
                        'Heterodimer dG','L Hairpin dG','R Hairpin dG',
                        'L Homodimer dG','R Homodimer dG','Thermo checks passed','Score']
        with open(str(csv_path), 'w', newline='') as _cf2:
            _cw2 = _csv_mod2.writer(_cf2)
            _cw2.writerow(_csv_hdrs_mc)
            _cw2.writerows(_csv_rows_mc)

        # Build subgenome groups
        from collections import defaultdict
        grouped = defaultdict(list)
        for _p in pairs:
            _sn = _p.seq_name.split()[0]; _p.seq_name = _sn; grouped[_sn].append(_p)

        from alignment import read_fasta as _rf3
        all_seq_names = [n.split()[0] for n in _rf3(str(fasta_path)).keys()]
        for _n in all_seq_names:
            if _n not in grouped:
                grouped[_n] = []

        subgenomes = []
        for _sn in all_seq_names:
            _sps = grouped.get(_sn, [])
            _sps.sort(key=lambda _p: (
                -(int(_p.left_genome_specific) + int(_p.right_genome_specific)),
                -(int(_p.left_end_genome_specific) + int(_p.right_end_genome_specific)),
                -_p.score))
            _top = []
            for _rk, _p in enumerate(_sps, 1):
                _top.append({
                    'left_genome_specific':  _p.left_genome_specific,
                    'right_genome_specific': _p.right_genome_specific,
                    'left_genome_hits':      _p.left_genome_hits,
                    'right_genome_hits':     _p.right_genome_hits,
                    'off_target_amplicon':   getattr(_p, 'off_target_amplicon', False),
                    'rank': _rk, 'seq_name': _p.seq_name,
                    'left_seq': _p.left_seq, 'left_start': _p.left_start,
                    'left_len': _p.left_len, 'left_tm': _p.left_tm, 'left_gc': _p.left_gc,
                    'right_seq': _p.right_seq, 'right_start': _p.right_start,
                    'right_len': _p.right_len, 'right_tm': _p.right_tm, 'right_gc': _p.right_gc,
                    'product_size': _p.product_size,
                    'left_diff': _p.left_diff_count, 'left_unique_pos': _p.left_unique_pos or [],
                    'right_diff': _p.right_diff_count, 'right_unique_pos': _p.right_unique_pos or [],
                    'left_3end': _p.left_end_genome_specific,
                    'right_3end': _p.right_end_genome_specific,
                    'score': round(_p.score, 2),
                    'thermo': {
                        'heterodimer':    calc_thermo(_p.left_seq, _p.right_seq, 'END1'),
                        'left_hairpin':   calc_thermo(_p.left_seq, mode='HAIRPIN'),
                        'right_hairpin':  calc_thermo(_p.right_seq, mode='HAIRPIN'),
                        'left_homodimer': calc_thermo(_p.left_seq, _p.left_seq, 'ANY'),
                        'right_homodimer':calc_thermo(_p.right_seq, _p.right_seq, 'ANY'),
                    },
                    'left_locations': [], 'right_locations': [],
                })
            subgenomes.append({'name': _sn, 'count': len(_sps), 'pairs': _top})

        msa_data = {
            'sequences':   [{'name': s.name.split()[0], 'aligned': s.aligned,
                             'length': len(s.aligned)} for s in msa.sequences],
            'unique_cols': list(set(u.col for u in unique)),
            'unique_by_seq': {s.name.split()[0]: sorted(u.col for u in unique if u.seq_name == s.name)
                              for s in msa.sequences},
            'msa_length':  msa.length,
        }
        top_pairs = [{'rank': _rk, 'seq_name': _p.seq_name,
                      'left_seq': _p.left_seq, 'left_start': _p.left_start,
                      'left_len': _p.left_len, 'left_tm': _p.left_tm, 'left_gc': _p.left_gc,
                      'right_seq': _p.right_seq, 'right_start': _p.right_start,
                      'right_len': _p.right_len, 'right_tm': _p.right_tm, 'right_gc': _p.right_gc,
                      'product_size': _p.product_size,
                      'left_diff': _p.left_diff_count, 'right_diff': _p.right_diff_count,
                      'score': round(_p.score, 2)}
                     for _rk, _p in enumerate(pairs, 1)]

        _set_job_done(job_id, {
            'success': True, 'db_name': db_name,
            'total_pairs': len(pairs),
            'msa_data': msa_data, 'top_pairs': top_pairs,
            'subgenomes': subgenomes, 'job_id': job_id,
            'warning': _repeat_warning,  # None if no repeat warning
        })

    except subprocess.TimeoutExpired:
        _set_job_error(job_id,
            'Alignment timed out — your sequence may be too long. '
            'Try submitting a shorter region (under 10,000 bp) around your SNP of interest.')
    except Exception as _exc:
        import traceback, logging as _logging
        _logging.getLogger(__name__).error('Design job %s failed: %s', job_id, traceback.format_exc())
        _set_job_error(job_id, 'An unexpected error occurred during primer design. Please try again or use a different sequence.')


# ── Progress polling endpoint ─────────────────────────────────────────────────

def job_progress(request, job_id: str):
    """Return current progress or final result for a background job."""
    _cleanup_old_jobs()
    job = _read_job(job_id)
    if job is None:
        return JsonResponse({'error': 'Job not found'}, status=404)

    elapsed = int(time.time() - job.get('started_at', time.time()))
    if job['status'] == 'done':
        result = dict(job['result'])
        result['elapsed'] = elapsed
        return JsonResponse(result)
    if job['status'] == 'error':
        return JsonResponse({'error': job['error'], 'elapsed': elapsed}, status=400)

    # Still running — return progress info
    return JsonResponse({
        'status':    'running',
        'step':      job.get('step', 1),
        'step_name': job.get('step_name', ''),
        'detail':    job.get('detail', ''),
        'elapsed':   elapsed,
        'steps':     job.get('steps_list', BLAST_STEPS),
    })


# ── Main view ─────────────────────────────────────────────────────────────────

@csrf_exempt
def run_pipeline(request):
    """Run BLAST + primer design pipeline."""
    if request.method != 'POST':
        return JsonResponse({'error': 'POST required'}, status=405)

    sequence   = request.POST.get('sequence', '').strip()
    fasta_file = request.FILES.get('fasta_file')
    mode       = request.POST.get('mode', 'blast')  # blast or direct
    db_name    = request.POST.get('database', '')
    db_path    = None

    # Initialise variables that are only set inside the BLAST branch so that
    # the code after the if/else can reference them safely without using dir().
    group     = None   # HitGroup — set only in BLAST mode
    extracted = {}     # {chrom: sequence} — set only in BLAST mode

    # Create a unique working directory
    job_id  = str(uuid.uuid4())[:8]
    workdir = Path(settings.MEDIA_ROOT) / 'jobs' / job_id
    workdir.mkdir(parents=True, exist_ok=True)

    # ── BLAST mode: validate input, start background thread, return immediately
    if mode == 'blast':
        if not sequence and not fasta_file:
            return JsonResponse({'error': 'No sequence provided'}, status=400)
        raw_seq_text = (fasta_file.read().decode('utf-8', errors='replace')
                        if fasta_file else sequence)
        if fasta_file:
            fasta_file.seek(0)
        err = _validate_fasta(raw_seq_text, min_length=100)
        if err:
            return JsonResponse({'error': err}, status=400)
        # BLAST mode only supports a single query sequence
        _n_seqs = sum(1 for _l in raw_seq_text.splitlines() if _l.strip().startswith('>'))
        if _n_seqs == 0:
            _n_seqs = 1  # bare sequence with no header
        if _n_seqs > 1:
            return JsonResponse(
                {'error': f'BLAST mode accepts one query sequence at a time '
                          f'(you provided {_n_seqs}). Please paste a single FASTA entry.'},
                status=400)
        db_path = settings.BLAST_DATABASES.get(db_name)
        if not db_path:
            return JsonResponse({'error': f'Unknown database: {db_name}'}, status=400)
        query_path = workdir / 'query.fas'
        if fasta_file:
            query_path.write_bytes(fasta_file.read())
        else:
            if not sequence.startswith('>'):
                sequence = f'>query\n{sequence}'
            query_path.write_text(sequence)
        amp_min_v = max(50,  min(int(request.POST.get('amp_min', 150)), 5000))
        amp_max_v = max(amp_min_v + 50, min(int(request.POST.get('amp_max', 700)), 5000))
        # Register job
        _write_job(job_id, {
            'status': 'running', 'step': 1,
            'step_name': BLAST_STEPS[0][1], 'detail': '',
            'result': None, 'error': None,
            'started_at': time.time(),
        })
        # Start background thread
        _t = _threading.Thread(
            target=_run_blast_pipeline,
            args=(job_id, workdir, query_path, db_path, db_name, amp_min_v, amp_max_v),
            daemon=True,
        )
        _t.start()
        return JsonResponse({
            'async':   True,
            'job_id':  job_id,
            'steps':   BLAST_STEPS,
        })

    # ── Direct mode: run synchronously (no genome BLAST, fast enough) ─────────
    try:
        if mode == 'direct':
            # Direct FASTA mode — supports multiple files
            fasta_files = request.FILES.getlist('fasta_files')
            if fasta_files:
                raw_bytes = b''
                for f in fasta_files:
                    raw_bytes += f.read() + b'\n'
                raw_text = raw_bytes.decode('utf-8', errors='replace')
                err = _validate_fasta(raw_text, min_length=100)
                if err:
                    return JsonResponse({'error': err}, status=400)
                fasta_path = workdir / 'input.fas'
                fasta_path.write_bytes(raw_bytes)
            elif fasta_file:
                raw_text = fasta_file.read().decode('utf-8', errors='replace')
                err = _validate_fasta(raw_text, min_length=100)
                if err:
                    return JsonResponse({'error': err}, status=400)
                fasta_path = workdir / 'input.fas'
                fasta_path.write_text(raw_text)
            elif sequence:
                err = _validate_fasta(sequence, min_length=100)
                if err:
                    return JsonResponse({'error': err}, status=400)
                fasta_path = workdir / 'input.fas'
                if not sequence.startswith('>'):
                    sequence = f'>query\n{sequence}'
                fasta_path.write_text(sequence)
            else:
                return JsonResponse({'error': 'No sequence provided'}, status=400)

        else:
            # Should not reach here — blast mode is handled above and returns early.
            return JsonResponse({'error': 'Unexpected mode'}, status=400)

        # Run primer design
        from alignment import build_msa
        from primers import design_primers, load_primer3_settings
        from output import write_csv

        # Require at least 2 sequences for the MSA + discriminating-position approach
        from alignment import read_fasta as _count_seqs
        _seq_count = len(_count_seqs(str(fasta_path)))
        if _seq_count < 2:
            return JsonResponse(
                {'error': 'Please provide at least 2 sequences so discriminating positions '
                          'can be identified. For a single sequence use the BLAST Against Genome tab.'},
                status=400)

        msa, unique = build_msa(None, str(fasta_path), settings.MUSCLE_PATH)

        if not unique:
            return JsonResponse({'error': 'No unique positions found — the sequences may be '
                                          'too similar to design discriminating primers. '
                                          'Check that you have provided the correct homeolog sequences.'}, status=400)

        # Apply user-specified amplicon size range
        amp_min = int(request.POST.get('amp_min', 150))
        amp_max = int(request.POST.get('amp_max', 1500))
        amp_min = max(50,  min(amp_min, 5000))
        amp_max = max(amp_min + 50, min(amp_max, 5000))
        from primers import DEFAULT_SETTINGS
        primer_settings = dict(DEFAULT_SETTINGS)
        primer_settings['PRIMER_PRODUCT_SIZE_RANGE'] = f'{amp_min}-{amp_max}'
        # Primer design parameters per specification
        primer_settings['PRIMER_MIN_GC']           = '35.0'
        primer_settings['PRIMER_MAX_GC']           = '65.0'
        primer_settings['PRIMER_MIN_TM']           = '50.0'
        primer_settings['PRIMER_OPT_TM']           = '60.0'
        primer_settings['PRIMER_MAX_TM']           = '65.0'
        primer_settings['PRIMER_PAIR_MAX_DIFF_TM'] = '2.0'
        primer_settings['PRIMER_MAX_POLY_X']       = '3'
        primer_settings['PRIMER_MIN_SIZE']         = '18'
        primer_settings['PRIMER_MAX_SIZE']         = '25'

        # FIX 1+2+3: pass full homeolog sequences and genome BLAST db
        # so primers are filtered genome-wide before pairing.
        from alignment import load_full_sequences_for_group

        # Build full_seqs for homeolog-level in-memory check
        # In both modes, fasta_path contains the sequences used for alignment.
        # read_fasta gives us the full sequences as written (no windowing here).
        from alignment import read_fasta as _read_fasta
        _full_seqs = _read_fasta(str(fasta_path))

        # target_chroms: ALL homeolog chromosomes selected by BLAST (e.g. 6A, 6B, 6D).
        # Hits on these are expected — the MSA uniqueness filter already handles
        # discrimination between homeologs. Only hits on completely different
        # chromosomes (1A, 1B, 2A ...) are true off-targets and disqualify a primer.
        # In direct mode there is no BLAST context so use None (strict 1-hit mode).
        _target_chroms = None
        if mode == 'blast' and group is not None:
            # Normalize chromosome names the same way primers.py normalizes sseqid:
            #   raw_id.split("|")[-1].split(":")[-1]
            # This ensures that Jagger-style IDs like "gnl|10Jagger|Ta_Jagger_7B"
            # become "Ta_Jagger_7B" and match what genome_blast_filter extracts
            # from BLAST output.  Without this, every primer fails the full-length
            # filter because "Ta_Jagger_7B" is never found in
            # {"gnl|10Jagger|Ta_Jagger_7B", ...}.
            _target_chroms = {
                n.split("|")[-1].split(":")[-1] for n in extracted.keys()
            }
            import logging as _lg
            _lg.getLogger(__name__).info("target_chroms (all homeologs): %s", _target_chroms)

        # query_chrom is the name of the target chromosome in the MSA.
        # MUSCLE reorders sequences so we cannot rely on position — pass it
        # explicitly via primary_target so design_primers finds the right sequence.
        _primary_target = query_chrom if mode == 'blast' else None

        pairs = design_primers(
            msa, unique, settings.PRIMER3_PATH,
            settings            = primer_settings,
            full_seqs           = _full_seqs,          # FIX 1: full homeolog check
            blast_db            = db_path if mode == 'blast' else None,  # NEW: genome-wide filter
            blastn_path         = settings.BLASTN_PATH,
            target_chroms       = _target_chroms,
            require_3end_unique = True,                # FIX 2: enforce 3' end
            primary_target      = _primary_target,     # FIX: MUSCLE reorders; pass explicit target name
        )

        if not pairs:
            return JsonResponse({'error': 'No specific primer pairs found'}, status=400)

        # For BLAST mode, only keep primers for the query sequence (first sequence)
        if mode == 'blast':
            query_seq_name = list(extracted.keys())[0] if extracted else None
            if query_seq_name:
                pairs = [p for p in pairs if p.seq_name == query_seq_name]
            if not pairs:
                return JsonResponse({'error': 'No specific primer pairs found for your query sequence'}, status=400)

        # Write multi-sheet Excel file
        import openpyxl
        from collections import defaultdict as _ddict
        xlsx_path = workdir / 'primers.xlsx'
        wb = openpyxl.Workbook()
        wb.remove(wb.active)  # remove default sheet

        headers = ['Rank','Sequence','Left primer (5->3)','Left start','Left len',
                   'Left Tm','Left GC%','Right primer (5->3)','Right start','Right len',
                   'Right Tm','Right GC%','Product size','Off-target amplicon',
                   'Heterodimer dG','L Hairpin dG','R Hairpin dG',
                   'L Homodimer dG','R Homodimer dG','Thermo checks passed','Score']

        # Group by seq_name
        by_seq = _ddict(list)
        for p in pairs:
            by_seq[p.seq_name].append(p)

        import re as _re2, csv as _csv_mod_t1
        _csv_rows_t1 = []
        for seq_name, seq_pairs in by_seq.items():
            chrom_match = _re2.search(r'([1-7][ABDSabds])', seq_name)
            if chrom_match:
                sheet_name = chrom_match.group(1).upper()
            else:
                sheet_name = _re2.sub(r'[:/\\?*\[\]]', '_', seq_name.split()[0])[:31]
            ws = wb.create_sheet(title=sheet_name)
            ws.append(headers)
            for p in seq_pairs:
                _th_t1 = {
                    'heterodimer':    calc_thermo(p.left_seq, p.right_seq, 'END1').get('dg'),
                    'left_hairpin':   calc_thermo(p.left_seq, mode='HAIRPIN').get('dg'),
                    'right_hairpin':  calc_thermo(p.right_seq, mode='HAIRPIN').get('dg'),
                    'left_homodimer': calc_thermo(p.left_seq, p.left_seq, 'ANY').get('dg'),
                    'right_homodimer':calc_thermo(p.right_seq, p.right_seq, 'ANY').get('dg'),
                }
                _rk_t1 = p.rank if hasattr(p, 'rank') else seq_pairs.index(p)+1
                _ot_t1 = 'Yes' if getattr(p, 'off_target_amplicon', False) else 'No'
                _row_t1 = [_rk_t1, p.seq_name, p.left_seq, p.left_start, p.left_len,
                            p.left_tm, p.left_gc, p.right_seq, p.right_start, p.right_len,
                            p.right_tm, p.right_gc, p.product_size, _ot_t1,
                            _thermo_risk_label(_th_t1['heterodimer']),
                            _thermo_risk_label(_th_t1['left_hairpin']),
                            _thermo_risk_label(_th_t1['right_hairpin']),
                            _thermo_risk_label(_th_t1['left_homodimer']),
                            _thermo_risk_label(_th_t1['right_homodimer']),
                            _thermo_summary([_th_t1['heterodimer'], _th_t1['left_hairpin'],
                                             _th_t1['right_hairpin'], _th_t1['left_homodimer'],
                                             _th_t1['right_homodimer']]),
                            round(p.score, 2)]
                ws.append(_row_t1)
                _csv_rows_t1.append(_row_t1)

        _autofit_xlsx(wb)
        wb.save(str(xlsx_path))

        # Write CSV in same format as xlsx
        csv_path = workdir / 'primers.csv'
        _csv_hdrs_t1 = ['Rank','Sequence','Left primer (5->3)','Left start','Left len',
                        'Left Tm','Left GC%','Right primer (5->3)','Right start','Right len',
                        'Right Tm','Right GC%','Product size','Off-target amplicon',
                        'Heterodimer dG','L Hairpin dG','R Hairpin dG',
                        'L Homodimer dG','R Homodimer dG','Thermo checks passed','Score']
        with open(str(csv_path), 'w', newline='') as _cf_t1:
            _cw_t1 = _csv_mod_t1.writer(_cf_t1)
            _cw_t1.writerow(_csv_hdrs_t1)
            _cw_t1.writerows(_csv_rows_t1)

        # Group pairs by subgenome (seq_name)
        from collections import defaultdict
        grouped = defaultdict(list)
        for p in pairs:
            # shorten name to first word only e.g. '2A dna:chromosome...' -> '2A'
            short_name = p.seq_name.split()[0]
            p.seq_name = short_name
            grouped[short_name].append(p)

        # Sorting is now done inside design_primers() using the full
        # specificity tier (genome-specific L+R > L-only > R-only > MSA 3'end > etc.)
        # No re-sort needed here.

        # Build per-subgenome top 20
        # Include ALL input sequences as tabs, even if 0 pairs
        from alignment import read_fasta
        # Use short names (first word) to match grouped keys
        all_seq_names = [name.split()[0] for name in read_fasta(str(fasta_path)).keys()]
        # Add sequences with 0 pairs to grouped
        for name in all_seq_names:
            if name not in grouped:
                grouped[name] = []

        subgenomes = []
        for seq_name in all_seq_names:  # preserve input order
            seq_pairs = grouped.get(seq_name, [])
            seq_pairs.sort(key=lambda p: (
                # Direct mode has no genome BLAST — sort by MSA 3'-end discriminating
                # positions first, then by total unique sites, then by score.
                -(int(p.left_3end_diff) + int(p.right_3end_diff)),
                -(p.left_diff_count + p.right_diff_count),
                -p.score
            ))
            top = []
            for rank, p in enumerate(seq_pairs, 1):
                top.append({
                    'left_genome_specific':  p.left_genome_specific,
                    'right_genome_specific': p.right_genome_specific,
                    'left_genome_hits':      p.left_genome_hits,
                    'right_genome_hits':     p.right_genome_hits,
                    'rank':         rank,  # rank is per-subgenome, starts from 1
                    'seq_name':     p.seq_name,
                    'left_seq':     p.left_seq,
                    'left_start':   p.left_start,
                    'left_len':     p.left_len,
                    'left_tm':      p.left_tm,
                    'left_gc':      p.left_gc,
                    'right_seq':    p.right_seq,
                    'right_start':  p.right_start,
                    'right_len':    p.right_len,
                    'right_tm':     p.right_tm,
                    'right_gc':     p.right_gc,
                    'product_size': p.product_size,
                    'left_diff':    p.left_diff_count,
                    'left_unique_pos': p.left_unique_pos or [],
                    'right_diff':   p.right_diff_count,
                    'right_unique_pos': p.right_unique_pos or [],
                    # In direct mode (no BLAST db) use the MSA 3'-end discriminating
                    # position flag — the genome-BLAST flag is always False here.
                    'left_3end':    p.left_3end_diff,
                    'right_3end':   p.right_3end_diff,
                    'score':        round(p.score, 2),
                    'thermo': {
                        'heterodimer':     calc_thermo(p.left_seq, p.right_seq, 'END1'),
                        'left_hairpin':    calc_thermo(p.left_seq, mode='HAIRPIN'),
                        'right_hairpin':   calc_thermo(p.right_seq, mode='HAIRPIN'),
                        'left_homodimer':  calc_thermo(p.left_seq, p.left_seq, 'ANY'),
                        'right_homodimer': calc_thermo(p.right_seq, p.right_seq, 'ANY'),
                    },
                    'left_locations':  [],
                    'right_locations': [],
                })
            subgenomes.append({
                'name':  seq_name,
                'count': len(seq_pairs),
                'pairs': top,
            })

        # Send MSA alignment data for visualization
        msa_data = {
            'sequences': [
                {
                    'name': s.name.split()[0],
                    'aligned': s.aligned,
                    'length': len(s.aligned),
                }
                for s in msa.sequences
            ],
            'unique_cols': list(set(u.col for u in unique)),
            'unique_by_seq': {
                s.name.split()[0]: sorted(
                    u.col for u in unique if u.seq_name == s.name
                )
                for s in msa.sequences
            },
            'msa_length': msa.length,
        }

        # Also keep flat top_pairs for visualization — send ALL pairs
        # Use enumerate() instead of pairs.index(p) to avoid O(n²) behaviour.
        top_pairs = []
        for rank, p in enumerate(pairs, 1):
            top_pairs.append({
                'rank':         rank,
                'seq_name':     p.seq_name,
                'left_seq':     p.left_seq,
                'left_start':   p.left_start,
                'left_len':     p.left_len,
                'left_tm':      p.left_tm,
                'left_gc':      p.left_gc,
                'right_seq':    p.right_seq,
                'right_start':  p.right_start,
                'right_len':    p.right_len,
                'right_tm':     p.right_tm,
                'right_gc':     p.right_gc,
                'product_size': p.product_size,
                'left_diff':    p.left_diff_count,
                'right_diff':   p.right_diff_count,
                'score':        round(p.score, 2),
            })

        return JsonResponse({
            'success':     True,
            'db_name':     db_name,
            'total_pairs': len(pairs),
            'msa_data':    msa_data,
            'top_pairs':   top_pairs,
            'subgenomes':  subgenomes,
            'job_id':      job_id,
        })

    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)


def detect_snps(request):
    """Detect SNPs between two sequences — fast, no primer design."""
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'gsp_python'))

    if request.method != 'POST':
        return JsonResponse({'error': 'POST required'}, status=405)

    seq1 = request.POST.get('seq1', '').strip()
    seq2 = request.POST.get('seq2', '').strip()

    _DNA_CHARS = set('ACGTNRYSWKMBDHVacgtnryswkmbdhv')

    def _looks_like_seq(text):
        """True if >70 % of chars are valid IUPAC nucleotide characters."""
        if not text:
            return False
        dna = sum(1 for c in text if c in _DNA_CHARS)
        return dna / len(text) > 0.70

    def split_multifasta(s):
        """Return list of (header, sequence) from a FASTA string.

        Handles:
        - Standard multi-line FASTA  (>name\\nseq...)
        - Inline-sequence FASTA      (>name ATGCATGC... all on one line)
        - Bare sequence with no header
        - Multi-FASTA (two records) pasted into a single box
        """
        records = []
        header = None
        bases = []
        for line in s.splitlines():
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(bases).upper()))
                rest = line[1:].strip()
                # Split "name" from any inline sequence on the same header line.
                # e.g. ">WTseq ATGCATGC..." — name="WTseq", inline_seq="ATGCATGC..."
                words = rest.split(None, 1)
                header = words[0] if words else f'seq{len(records)+1}'
                bases = []
                if len(words) > 1:
                    candidate = words[1].replace(' ', '')
                    if _looks_like_seq(candidate):
                        bases.append(candidate)   # inline sequence on same line
                    # else: it's a description/annotation — ignore
            elif line:
                bases.append(line)
        if header is not None:
            records.append((header, ''.join(bases).upper()))
        # No FASTA header at all — treat entire string as a bare sequence
        if not records and s:
            bare = s.replace('\n', '').replace(' ', '').upper()
            records.append(('seq1', bare))
        return records

    def strip_fasta(s):
        """Extract first sequence from a FASTA or bare-sequence string."""
        recs = split_multifasta(s)
        return recs[0][1] if recs else ''

    # If seq2 is empty but seq1 contains two FASTA records, auto-split
    if not seq2:
        recs1 = split_multifasta(seq1)
        if len(recs1) >= 2:
            seq1 = recs1[0][1]
            seq2 = recs1[1][1]
        else:
            seq1 = recs1[0][1] if recs1 else ''
    else:
        seq1 = strip_fasta(seq1)
        seq2 = strip_fasta(seq2)

    if not seq1 or not seq2:
        return JsonResponse({'error': 'Both sequences required'}, status=400)
    if len(seq1) < 50 or len(seq2) < 50:
        return JsonResponse({'error': 'Sequences too short (min 50bp)'}, status=400)
    if seq1.upper() == seq2.upper():
        return JsonResponse({'error': 'The two sequences are identical — no SNPs can be detected.'}, status=400)

    try:
        from kasp import find_snps, align_two_seqs

        if len(seq1) != len(seq2):
            seq1, seq2 = align_two_seqs(seq1, seq2, settings.MUSCLE_PATH)

        snps = find_snps(seq1, seq2)
        return JsonResponse({
            # Report 1-based positions (biologists expect position 1 = first base)
            'snps': [{'pos': s.raw_pos + 1, 'allele1': s.allele1, 'allele2': s.allele2} for s in snps],
            'total': len(snps),
        })
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)

def kasp_design(request):
    """KASP marker design endpoint."""
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent.parent / 'gsp_python'))

    if request.method != 'POST':
        return JsonResponse({'error': 'POST required'}, status=405)

    seq1 = request.POST.get('seq1', '').strip()
    seq2 = request.POST.get('seq2', '').strip()
    snp_positions_raw = request.POST.get('snp_positions', '')
    database = request.POST.get('database', '').strip()
    # database is optional

    def split_multifasta(s):
        """Return list of (header, sequence) tuples from a FASTA or bare-sequence string."""
        _DNA_CHARS = set('ACGTNRYSWKMBDHVacgtnryswkmbdhv')

        def _looks_like_seq(text):
            if not text:
                return False
            dna = sum(1 for c in text if c in _DNA_CHARS)
            return dna / len(text) > 0.70

        records = []
        header = None
        bases = []
        for line in s.splitlines():
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(bases).upper()))
                rest = line[1:].strip()
                words = rest.split(None, 1)
                header = words[0] if words else f'seq{len(records)+1}'
                bases = []
                if len(words) > 1:
                    candidate = words[1].replace(' ', '')
                    if _looks_like_seq(candidate):
                        bases.append(candidate)
            elif line:
                bases.append(line)
        if header is not None:
            records.append((header, ''.join(bases).upper()))
        if not records and s:
            bare = s.replace('\n', '').replace(' ', '').upper()
            records.append(('seq1', bare))
        return records

    def strip_fasta(s):
        recs = split_multifasta(s)
        return recs[0][1] if recs else ''

    # If seq2 is empty but seq1 contains two FASTA records, auto-split
    if not seq2:
        recs1 = split_multifasta(seq1)
        if len(recs1) >= 2:
            seq1 = recs1[0][1]
            seq2 = recs1[1][1]
        else:
            seq1 = recs1[0][1] if recs1 else ''
    else:
        seq1 = strip_fasta(seq1)
        seq2 = strip_fasta(seq2)

    if not seq1 or not seq2:
        return JsonResponse({'error': 'Both sequences required — paste one sequence per box, or paste both as a multi-FASTA in the Parent 1 box.'}, status=400)

    if len(seq1) < 150 or len(seq2) < 150:
        return JsonResponse({'error': f'Sequences too short (min 150bp). Got {len(seq1)}bp and {len(seq2)}bp.'}, status=400)

    if seq1.upper() == seq2.upper():
        return JsonResponse({'error': 'Parent 1 and Parent 2 sequences are identical — no SNPs can be detected. Please check that you have pasted different sequences.'}, status=400)

    snp_positions = None
    if snp_positions_raw.strip():
        try:
            # Frontend sends 1-based positions (as displayed to the user).
            # kasp.py's find_snps uses 0-based raw_pos internally, so subtract 1.
            snp_positions = [int(x.strip()) - 1 for x in snp_positions_raw.split(',') if x.strip()]
        except:
            pass

    db_path = settings.BLAST_DATABASES.get(database, '')

    # ── Inner BLAST helpers (closed over db_path / settings) ─────────────────
    def check_specificity(primer_seq, get_locations=False):
        if not db_path:
            return (None, []) if get_locations else None
        try:
            import subprocess as _sp
            r = _sp.run([
                settings.BLASTN_PATH, '-db', db_path,
                '-query', '-', '-task', 'blastn-short',
                '-word_size', '11', '-outfmt', '6',
                '-max_target_seqs', '100', '-evalue', '1000',
                '-dust', 'no',
            ], input=f'>p\n{primer_seq}\n',
            capture_output=True, text=True, timeout=30)
            chroms = set()
            locations = []
            seen_chroms = set()
            for hit in r.stdout.strip().split('\n'):
                if hit:
                    parts = hit.split('\t')
                    if len(parts) >= 9 and float(parts[2]) >= 95 and int(parts[3]) >= max(15, int(0.85 * len(primer_seq))):
                        chrom = parts[1]
                        chrom_short = chrom.split(':')[0]
                        chroms.add(chrom_short)
                        if chrom_short not in seen_chroms:
                            seen_chroms.add(chrom_short)
                            locations.append({
                                'chrom': chrom_short,
                                'start': int(parts[8]),
                                'end': int(parts[9]),
                                'identity': float(parts[2]),
                            })
            if get_locations:
                return len(chroms), locations
            return len(chroms)
        except Exception:
            return (None, []) if get_locations else None

    def get_blast_hits_with_coords(primer_seq):
        """Return list of (chrom, sstart, send, sstrand, pident, aln_len) for significant hits."""
        if not db_path:
            return []
        try:
            import subprocess as _sp
            r = _sp.run([
                settings.BLASTN_PATH, '-db', db_path,
                '-query', '-', '-task', 'blastn-short',
                '-word_size', '7',
                '-outfmt', '6 qseqid sseqid pident length qlen qstart qend sstart send sstrand',
                '-max_target_seqs', '200', '-evalue', '1',
                '-dust', 'no', '-perc_identity', '80',
            ], input=f'>p\n{primer_seq}\n',
            capture_output=True, text=True, timeout=30)
            hits = []
            n = len(primer_seq)
            for line in r.stdout.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 10:
                    continue
                try:
                    pident  = float(parts[2])
                    aln_len = int(parts[3])
                    sseqid  = parts[1]
                    sstart  = int(parts[7])
                    send    = int(parts[8])
                    sstrand = parts[9].strip()
                except (ValueError, IndexError):
                    continue
                if pident >= 95.0 and aln_len >= max(15, int(0.85 * n)):
                    chrom = sseqid.split(':')[0]
                    hits.append((chrom, sstart, send, sstrand, pident, aln_len))
            return hits
        except Exception:
            return []

    def check_pair_offtarget(fwd_primer, rev_primer, max_amplicon=170):
        """Return True if fwd+rev can form an off-target amplicon (50–max_amplicon bp)."""
        fwd_hits = get_blast_hits_with_coords(fwd_primer)
        rev_hits = get_blast_hits_with_coords(rev_primer)
        if not fwd_hits or not rev_hits:
            return False
        fwd_by_chr: dict = {}
        for chrom, sstart, send, sstrand, *_ in fwd_hits:
            fwd_by_chr.setdefault(chrom, []).append((sstart, send, sstrand))
        rev_by_chr: dict = {}
        for chrom, sstart, send, sstrand, *_ in rev_hits:
            rev_by_chr.setdefault(chrom, []).append((sstart, send, sstrand))
        for chrom in fwd_by_chr:
            if chrom not in rev_by_chr:
                continue
            for f_start, f_end, f_strand in fwd_by_chr[chrom]:
                for r_start, r_end, r_strand in rev_by_chr[chrom]:
                    if f_strand == 'plus' and r_strand == 'minus':
                        amplicon = r_end - f_start + 1
                    elif f_strand == 'minus' and r_strand == 'plus':
                        amplicon = f_end - r_start + 1
                    else:
                        continue
                    if 50 <= amplicon <= max_amplicon:
                        return True
        return False

    # ── Create async job and start background thread ───────────────────────────
    job_id = str(uuid.uuid4())[:8]
    _write_job(job_id, {
        'status':     'running',
        'step':       1,
        'step_name':  KASP_STEPS[0][1],
        'detail':     '',
        'started_at': time.time(),
        'steps_list': KASP_STEPS,
    })

    def _worker():
        try:
            # Step 1: Align sequences
            _kasp_progress(job_id, 1)
            from kasp import design_kasp_markers, FAM_TAIL, HEX_TAIL, align_two_seqs
            aln_seq1, aln_seq2 = align_two_seqs(seq1, seq2, settings.MUSCLE_PATH)

            # Step 2: Detect SNPs
            _kasp_progress(job_id, 2)

            # Step 3: Design primer candidates (Primer3, per-SNP/orientation updates via callback)
            _kasp_progress(job_id, 3)
            all_snps, markers = design_kasp_markers(
                aln_seq1, aln_seq2,
                primer3_bin=settings.PRIMER3_PATH,
                muscle_bin=settings.MUSCLE_PATH,
                snp_positions=snp_positions,
                max_markers=50,
                db_path=db_path,
                blastn_path=settings.BLASTN_PATH,
                blastdbcmd_path=settings.BLASTDBCMD_PATH,
                progress_callback=lambda detail: _kasp_detail(job_id, detail),
            )

            # Step 4: Genome BLAST specificity check
            _kasp_progress(job_id, 4)
            for m in markers:
                m.body_hits, m.body_locations = check_specificity(m.fwd1_seq, get_locations=True)
                m.rev_hits,  m.rev_locations  = check_specificity(m.rev_seq,  get_locations=True)

            # Step 5: Off-target amplicon check & thermodynamics
            _kasp_progress(job_id, 5)
            for m in markers:
                m.fam_offtarget_amplicon = check_pair_offtarget(m.fwd1_seq, m.rev_seq)
                m.hex_offtarget_amplicon = check_pair_offtarget(m.fwd2_seq, m.rev_seq)
                m.thermo = {
                    'fam_common':    calc_thermo(m.fwd1_seq, m.rev_seq, 'END1'),
                    'hex_common':    calc_thermo(m.fwd2_seq, m.rev_seq, 'END1'),
                    'fam_hairpin':   calc_thermo(m.fwd1_seq, mode='HAIRPIN'),
                    'hex_hairpin':   calc_thermo(m.fwd2_seq, mode='HAIRPIN'),
                    'rev_hairpin':   calc_thermo(m.rev_seq,  mode='HAIRPIN'),
                    'fam_homodimer': calc_thermo(m.fwd1_seq, m.fwd1_seq, 'ANY'),
                    'hex_homodimer': calc_thermo(m.fwd2_seq, m.fwd2_seq, 'ANY'),
                    'rev_homodimer': calc_thermo(m.rev_seq,  m.rev_seq,  'ANY'),
                }

            # Step 6: Finalise — sort, build Excel, assemble result
            _kasp_progress(job_id, 6)

            # Sort markers by SNP position then specificity tier
            if db_path:
                def _marker_sort_key(m):
                    rh = m.rev_hits  if m.rev_hits  is not None else 9999
                    bh = m.body_hits if m.body_hits is not None else 9999
                    if   rh <= 1: tier = 0
                    elif rh <= 3: tier = 1
                    else:         tier = 2
                    return (m.snp_pos, tier, rh, bh)
                markers.sort(key=_marker_sort_key)

            # Build Excel download
            _kasp_workdir = Path(settings.MEDIA_ROOT) / 'jobs' / job_id
            _kasp_workdir.mkdir(parents=True, exist_ok=True)
            import openpyxl as _opx
            _wb = _opx.Workbook(); _ws = _wb.active; _ws.title = 'KASP Markers'
            _ws.append(['SNP pos','Allele1','Allele2','Orientation','Product (bp)',
                'FAM gene-specific','FAM full (with tail)',
                'HEX gene-specific','HEX full (with tail)',
                'Common primer','FAM Tm','HEX Tm','Common Tm',
                'FAM GC%','HEX GC%','Common GC%'])
            for _m in markers:
                _ws.append([_m.snp_pos + 1,_m.allele1,_m.allele2,_m.orientation,_m.product_size,
                    _m.fwd1_seq,_m.fwd1_with_tail,_m.fwd2_seq,_m.fwd2_with_tail,_m.rev_seq,
                    round(_m.fwd1_tm,1),round(_m.fwd2_tm,1),round(_m.rev_tm,1),
                    round(_m.fwd1_gc,1),round(_m.fwd2_gc,1),round(_m.rev_gc,1)])
            _autofit_xlsx(_wb)
            _wb.save(str(_kasp_workdir / 'kasp_markers.xlsx'))

            # MSA viewer data
            msa_data = {
                'sequences': [
                    {'name': 'Parent1', 'aligned': aln_seq1, 'length': len(aln_seq1)},
                    {'name': 'Parent2', 'aligned': aln_seq2, 'length': len(aln_seq2)},
                ],
                'msa_length': len(aln_seq1),
                'unique_cols': [s.position for s in all_snps],
                'unique_by_seq': {
                    'Parent1': [s.position for s in all_snps],
                    'Parent2': [s.position for s in all_snps],
                },
            }

            result = {
                'job_id': job_id,
                'snps': [{'pos': s.raw_pos + 1, 'allele1': s.allele1, 'allele2': s.allele2} for s in all_snps],
                'markers': [{
                    'snp_pos':             m.snp_pos + 1,
                    'allele1':             m.allele1,
                    'allele2':             m.allele2,
                    'fwd1_seq':            m.fwd1_seq,
                    'fwd2_seq':            m.fwd2_seq,
                    'rev_seq':             m.rev_seq,
                    'fwd1_with_tail':      m.fwd1_with_tail,
                    'fwd2_with_tail':      m.fwd2_with_tail,
                    'fwd1_tm':             round(m.fwd1_tm, 1),
                    'fwd2_tm':             round(m.fwd2_tm, 1),
                    'rev_tm':              round(m.rev_tm, 1),
                    'fwd1_gc':             round(m.fwd1_gc, 1),
                    'fwd2_gc':             round(m.fwd2_gc, 1),
                    'rev_gc':              round(m.rev_gc, 1),
                    'product_size':        m.product_size,
                    'orientation':         m.orientation,
                    'fam_tail':            FAM_TAIL,
                    'hex_tail':            HEX_TAIL,
                    'body_hits':           getattr(m, 'body_hits', None),
                    'body_locations':      getattr(m, 'body_locations', []),
                    'rev_hits':            getattr(m, 'rev_hits', None),
                    'rev_locations':       getattr(m, 'rev_locations', []),
                    'fam_offtarget_amplicon': getattr(m, 'fam_offtarget_amplicon', None),
                    'hex_offtarget_amplicon': getattr(m, 'hex_offtarget_amplicon', None),
                    'thermo':              getattr(m, 'thermo', {}),
                } for m in markers],
                'total_snps':    len(all_snps),
                'total_markers': len(markers),
                'msa_data':      msa_data,
            }

            _update_job(job_id, status='done', result=result)

        except subprocess.TimeoutExpired:
            _update_job(job_id, status='error',
                        error='Design timed out — try fewer SNP positions or a smaller genome database')
        except Exception as e:
            import traceback as _tb
            _update_job(job_id, status='error', error=str(e), trace=_tb.format_exc())

    _threading.Thread(target=_worker, daemon=True).start()
    return JsonResponse({'async': True, 'job_id': job_id, 'steps': KASP_STEPS})


def download_selected(request):
    """Download selected primers as multi-sheet Excel."""
    import json, openpyxl, re
    from io import BytesIO

    if request.method != 'POST':
        return HttpResponse('Method not allowed', status=405)

    try:
        data = json.loads(request.body)
        pairs = data.get('pairs', [])
    except:
        return HttpResponse('Invalid data', status=400)

    if not pairs:
        return HttpResponse('No pairs', status=400)

    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    headers = ['Rank','Sequence','Left primer (5->3)','Left start','Left len',
               'Left Tm','Left GC%','Right primer (5->3)','Right start','Right len',
               'Right Tm','Right GC%','Product size','Off-target amplicon',
               'Heterodimer dG','L Hairpin dG','R Hairpin dG',
               'L Homodimer dG','R Homodimer dG','Thermo checks passed','Score']

    # Group by seq_name
    from collections import defaultdict
    grouped = defaultdict(list)
    for p in pairs:
        grouped[p['seq_name']].append(p)

    for seq_name, seq_pairs in grouped.items():
        chrom_match = re.search(r'([1-7][ABDS])\b', seq_name)
        sheet_name = chrom_match.group(1) if chrom_match else seq_name.split()[0][:31]
        ws = wb.create_sheet(title=sheet_name)
        ws.append(headers)
        for p in seq_pairs:
            th = p.get('thermo') or {}
            def _dg(key): return (th.get(key) or {}).get('dg')
            _ot_sel = 'Yes' if p.get('off_target_amplicon') else 'No'
            _dg_vals = [_dg('heterodimer'), _dg('left_hairpin'), _dg('right_hairpin'),
                        _dg('left_homodimer'), _dg('right_homodimer')]
            ws.append([
                p.get('rank'), p.get('seq_name'),
                p.get('left_seq'), p.get('left_start'), p.get('left_len'),
                p.get('left_tm'), p.get('left_gc'),
                p.get('right_seq'), p.get('right_start'), p.get('right_len'),
                p.get('right_tm'), p.get('right_gc'),
                p.get('product_size'), _ot_sel,
                _thermo_risk_label(_dg_vals[0]), _thermo_risk_label(_dg_vals[1]),
                _thermo_risk_label(_dg_vals[2]), _thermo_risk_label(_dg_vals[3]),
                _thermo_risk_label(_dg_vals[4]),
                _thermo_summary(_dg_vals),
                p.get('score')
            ])

    buf = BytesIO()
    _autofit_xlsx(wb)
    wb.save(buf)
    buf.seek(0)

    response = HttpResponse(buf.read(),
        content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    response['Content-Disposition'] = 'attachment; filename="selected_primers.xlsx"'
    return response

def download_csv(request, job_id):
    """Download primers as Excel or CSV based on fmt param."""
    fmt = request.GET.get('fmt', 'xlsx')
    xlsx_path = Path(settings.MEDIA_ROOT) / 'jobs' / job_id / 'primers.xlsx'
    csv_path  = Path(settings.MEDIA_ROOT) / 'jobs' / job_id / 'primers.csv'
    if fmt == 'csv' and csv_path.exists():
        response = HttpResponse(csv_path.read_bytes(), content_type='text/csv')
        response['Content-Disposition'] = f'attachment; filename="primers_{job_id}.csv"'
        return response
    if xlsx_path.exists():
        response = HttpResponse(xlsx_path.read_bytes(),
            content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = f'attachment; filename="primers_{job_id}.xlsx"'
        return response
    elif csv_path.exists():
        response = HttpResponse(csv_path.read_bytes(), content_type='text/csv')
        response['Content-Disposition'] = f'attachment; filename="primers_{job_id}.csv"'
        return response
    return HttpResponse('File not found', status=404)


def kasp_download_all(request, job_id):
    from pathlib import Path
    xlsx = Path(settings.MEDIA_ROOT) / 'jobs' / job_id / 'kasp_markers.xlsx'
    if not xlsx.exists():
        return HttpResponse('File not found', status=404)
    response = HttpResponse(xlsx.read_bytes(),
        content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    response['Content-Disposition'] = f'attachment; filename="kasp_markers_{job_id}.xlsx"'
    return response


def kasp_download_selected(request):
    import json, openpyxl
    from io import BytesIO
    if request.method != 'POST':
        return HttpResponse('Method not allowed', status=405)
    try:
        data = json.loads(request.body)
        markers = data.get('markers', [])
    except:
        return HttpResponse('Invalid data', status=400)
    if not markers:
        return HttpResponse('No markers', status=400)
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = 'Selected KASP Markers'
    ws.append(['SNP pos','Allele1','Allele2','Orientation','Product (bp)',
        'FAM gene-specific','FAM full (with tail)',
        'HEX gene-specific','HEX full (with tail)',
        'Common primer','FAM Tm','HEX Tm','Common Tm',
        'FAM GC%','HEX GC%','Common GC%'])
    for m in markers:
        ws.append([m.get('snp_pos'),m.get('allele1'),m.get('allele2'),
            m.get('orientation'),m.get('product_size'),
            m.get('fwd1_seq'),m.get('fwd1_with_tail'),
            m.get('fwd2_seq'),m.get('fwd2_with_tail'),m.get('rev_seq'),
            m.get('fwd1_tm'),m.get('fwd2_tm'),m.get('rev_tm'),
            m.get('fwd1_gc'),m.get('fwd2_gc'),m.get('rev_gc')])
    _autofit_xlsx(wb)
    buf = BytesIO(); wb.save(buf); buf.seek(0)
    response = HttpResponse(buf.read(),
        content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    response['Content-Disposition'] = 'attachment; filename="selected_kasp_markers.xlsx"'
    return response


@csrf_exempt
def primer_locations(request):
    """Fetch chromosomal locations for a single primer pair on demand."""
    import json
    if request.method != 'POST':
        return JsonResponse({'error': 'POST required'}, status=405)
    try:
        data       = json.loads(request.body)
        left_seq   = data.get('left_seq', '')
        right_seq  = data.get('right_seq', '')
        db_name    = data.get('db_name', '')
    except:
        return JsonResponse({'error': 'Invalid data'}, status=400)

    db_path = settings.BLAST_DATABASES.get(db_name, '')
    if not db_path:
        return JsonResponse({'error': f'Unknown database: {db_name}'}, status=400)
    # Run each BLAST once and reuse the results (was 6 calls, now 3).
    # left_full_specific / right_full_specific: exactly 1 hit = specific (was wrongly >= 1).
    left_locs  = _primer_locations(left_seq,  db_path, max_locs=10)
    right_locs = _primer_locations(right_seq, db_path, max_locs=10)
    left_3end  = _primer_3end_locations(left_seq,  db_path)
    right_3end = _primer_3end_locations(right_seq, db_path)
    return JsonResponse({
        'left_locations':      left_locs,
        'right_locations':     right_locs,
        'left_3end_locs':      left_3end,
        'right_3end_locs':     right_3end,
        'left_3end_specific':  len(left_3end)  == 1,
        'right_3end_specific': len(right_3end) == 1,
        'left_full_specific':  len(left_locs)  == 1,
        'right_full_specific': len(right_locs) == 1,
        'pair_safe':           _check_pair_off_target(left_seq, right_seq, db_path),
    })
