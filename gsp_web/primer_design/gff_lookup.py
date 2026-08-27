"""
gff_lookup.py
Tabix-based queries against bgzipped, CSI-indexed GFF3 files.

Public API
----------
lookup_gene_by_id(gene_id, gff_paths)        -> dict | None
search_genes_by_keyword(keyword, chrom, gff_paths, max_results=50) -> list[dict]
search_genes_by_region(chrom, start, end, gff_paths, max_results=50) -> list[dict]
get_cds_intervals(gene_id, chrom, gff_paths) -> list[(start, end, strand)]
get_intron_intervals(gene_id, chrom, gff_paths) -> list[(start, end, strand)]
extract_sequence(chrom, start, end, strand, blastdbcmd_path, blast_db) -> str
"""

from __future__ import annotations
import re
import subprocess
import logging
from typing import Optional

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_attrs(attr_str: str) -> dict[str, str]:
    """Parse GFF3 attribute column into a dict."""
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(';'):
        part = part.strip()
        if '=' in part:
            k, v = part.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs


def _gff_row_to_dict(cols: list[str]) -> dict:
    attrs = _parse_attrs(cols[8])
    gene_id = attrs.get('ID', '')
    # strip transcript suffix if present (e.g. TraesJAG1A01G000100.1 → TraesJAG1A01G000100)
    base_id = gene_id.split('.')[0] if '.' in gene_id else gene_id
    note = attrs.get('Note', '')
    name = attrs.get('Name', base_id)
    return {
        'gene_id':   base_id,
        'chrom':     cols[0],
        'start':     int(cols[3]),
        'end':       int(cols[4]),
        'strand':    cols[6],
        'feature':   cols[2],
        'note':      note,
        'name':      name,
        'attrs':     attrs,
    }


def _open_tabix(gff_path: str):
    """Return a pysam.TabixFile, or None if unavailable.
    Tries CSI index first (used by wheat genome GFFs), falls back to TBI."""
    try:
        import pysam  # type: ignore
        csi_path = gff_path + '.csi'
        tbi_path = gff_path + '.tbi'
        import os
        if os.path.exists(csi_path):
            return pysam.TabixFile(gff_path, index=csi_path)
        elif os.path.exists(tbi_path):
            return pysam.TabixFile(gff_path)
        else:
            log.warning("No tabix index (.csi or .tbi) found for %s", gff_path)
            return None
    except Exception as e:
        log.warning("Cannot open tabix file %s: %s", gff_path, e)
        return None


def _grep_gene_id(gene_id: str, gff_path: str) -> Optional[dict]:
    """
    Scan a bgzipped GFF3 to find a gene by ID.
    Uses Python's gzip module — works on macOS and Linux.
    """
    import gzip
    try:
        with gzip.open(gff_path, 'rt', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9 or cols[2] != 'gene':
                    continue
                attrs = _parse_attrs(cols[8])
                gid = attrs.get('ID', '').split('.')[0]
                if gid == gene_id:
                    return _gff_row_to_dict(cols)
    except Exception as e:
        log.warning("gzip scan failed for %s: %s", gene_id, e)
    return None


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

def lookup_gene_by_id(gene_id: str, gff_paths: list[str]) -> Optional[dict]:
    """
    Find a gene by its ID across a list of GFF3 files.
    Returns a dict with gene info, or None if not found.
    """
    gene_id = gene_id.strip()
    for gff_path in gff_paths:
        result = _grep_gene_id(gene_id, gff_path)
        if result:
            return result
    return None


def search_genes_by_keyword(
    keyword: str,
    chrom: str,
    gff_paths: list[str],
    max_results: int = 50,
) -> list[dict]:
    """
    Search genes on a chromosome whose Note/Name attributes contain keyword.
    Uses tabix for fast chromosome-level fetch, then filters by keyword.
    """
    keyword_lower = keyword.lower()
    results: list[dict] = []

    for gff_path in gff_paths:
        tbx = _open_tabix(gff_path)
        if tbx is None:
            continue
        try:
            for row in tbx.fetch(chrom):
                cols = row.split('\t')
                if len(cols) < 9 or cols[2] != 'gene':
                    continue
                attrs = _parse_attrs(cols[8])
                searchable = (
                    attrs.get('Note', '') + ' ' +
                    attrs.get('Name', '') + ' ' +
                    attrs.get('ID', '')
                ).lower()
                if keyword_lower in searchable:
                    results.append(_gff_row_to_dict(cols))
                    if len(results) >= max_results:
                        return results
        except Exception as e:
            log.warning("tabix fetch failed for %s on %s: %s", gff_path, chrom, e)
        finally:
            tbx.close()

    return results


def search_genes_by_region(
    chrom: str,
    start: int,
    end: int,
    gff_paths: list[str],
    max_results: int = 50,
) -> list[dict]:
    """
    Return all genes in a genomic region (1-based, inclusive).
    """
    results: list[dict] = []

    for gff_path in gff_paths:
        tbx = _open_tabix(gff_path)
        if tbx is None:
            continue
        try:
            for row in tbx.fetch(chrom, start - 1, end):  # pysam is 0-based
                cols = row.split('\t')
                if len(cols) < 9 or cols[2] != 'gene':
                    continue
                results.append(_gff_row_to_dict(cols))
                if len(results) >= max_results:
                    return results
        except Exception as e:
            log.warning("tabix region fetch failed: %s", e)
        finally:
            tbx.close()

    return results


def get_cds_intervals(
    gene_id: str,
    chrom: str,
    gff_paths: list[str],
) -> list[tuple[int, int, str]]:
    """
    Return sorted CDS intervals for a gene as (start, end, strand) tuples (1-based).
    Intervals are sorted by start position; for minus-strand genes, caller should
    reverse-complement the concatenated sequence.
    """
    intervals: list[tuple[int, int, str]] = []

    for gff_path in gff_paths:
        tbx = _open_tabix(gff_path)
        if tbx is None:
            continue
        try:
            for row in tbx.fetch(chrom):
                cols = row.split('\t')
                if len(cols) < 9 or cols[2] != 'CDS':
                    continue
                attrs = _parse_attrs(cols[8])
                parent = attrs.get('Parent', '').split('.')[0]
                if parent == gene_id:
                    intervals.append((int(cols[3]), int(cols[4]), cols[6]))
        except Exception as e:
            log.warning("CDS fetch failed: %s", e)
        finally:
            tbx.close()
        if intervals:
            break  # found in this file

    return sorted(intervals, key=lambda x: x[0])


def get_intron_intervals(
    gene_id: str,
    chrom: str,
    gff_paths: list[str],
) -> list[tuple[int, int, str]]:
    """
    Derive intron intervals from the gaps between sorted exons.
    Returns (start, end, strand) tuples (1-based).
    """
    exons: list[tuple[int, int, str]] = []

    for gff_path in gff_paths:
        tbx = _open_tabix(gff_path)
        if tbx is None:
            continue
        try:
            for row in tbx.fetch(chrom):
                cols = row.split('\t')
                if len(cols) < 9 or cols[2] != 'exon':
                    continue
                attrs = _parse_attrs(cols[8])
                parent = attrs.get('Parent', '').split('.')[0]
                if parent == gene_id:
                    exons.append((int(cols[3]), int(cols[4]), cols[6]))
        except Exception as e:
            log.warning("Exon fetch failed: %s", e)
        finally:
            tbx.close()
        if exons:
            break

    exons = sorted(exons, key=lambda x: x[0])
    introns: list[tuple[int, int, str]] = []
    for i in range(len(exons) - 1):
        intron_start = exons[i][1] + 1
        intron_end   = exons[i + 1][0] - 1
        strand       = exons[i][2]
        if intron_start <= intron_end:
            introns.append((intron_start, intron_end, strand))

    return introns


def extract_sequence(
    chrom: str,
    start: int,
    end: int,
    strand: str,
    blastdbcmd_path: str,
    blast_db: str,
) -> str:
    """
    Extract a sequence from a BLAST database using blastdbcmd.
    Coordinates are 1-based inclusive.
    Always extracts the plus-strand sequence; caller handles orientation.
    Returns the sequence string (uppercase), or '' on failure.
    """
    cmd = [
        blastdbcmd_path,
        '-db', blast_db,
        '-entry', chrom,
        '-range', f'{start}-{end}',
        '-strand', 'plus',
        '-outfmt', '%s',
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if result.returncode != 0:
            log.error("blastdbcmd error (rc=%d) for %s:%d-%d: %s",
                      result.returncode, chrom, start, end, result.stderr.strip())
            return ''
        seq = result.stdout.strip().upper()
        # Strip any accidental FASTA header lines
        if seq.startswith('>'):
            seq = ''.join(seq.split('\n')[1:])
        return seq
    except Exception as e:
        log.error("blastdbcmd failed for %s:%d-%d: %s", chrom, start, end, e)
        return ''


def _revcomp(seq: str) -> str:
    comp = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(comp)[::-1]


def list_chromosomes(gff_paths: list[str]) -> list[str]:
    """Return unique chromosome/seqid names present in the GFF file(s)."""
    chroms: set[str] = set()
    for gff_path in gff_paths:
        tbx = _open_tabix(gff_path)
        if tbx is None:
            continue
        try:
            chroms.update(tbx.contigs)
        except Exception:
            pass
        finally:
            tbx.close()
    return sorted(chroms)
