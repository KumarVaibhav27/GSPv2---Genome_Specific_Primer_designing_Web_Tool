"""
blast.py — Parse BLAST tabular output (format 6) into hit groups.
No changes from original — included for completeness.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class BlastHit:
    """One line from a BLAST tab file (format 6, 12 columns)."""
    qseqid:   str
    sseqid:   str
    pident:   float
    length:   int
    mismatch: int
    gapopen:  int
    qstart:   int
    qend:     int
    sstart:   int
    send:     int
    evalue:   float
    bitscore: float

    @property
    def strand(self) -> str:
        return "+" if self.sstart <= self.send else "-"

    @property
    def chrom_start(self) -> int:
        """0-based start for BED."""
        return min(self.sstart, self.send) - 1

    @property
    def chrom_end(self) -> int:
        """1-based end for BED."""
        return max(self.sstart, self.send)


@dataclass
class HitGroup:
    """All BLAST hits for one query sequence, grouped by subject chromosome."""
    query_id: str
    hits_by_chrom: dict[str, list[BlastHit]] = field(default_factory=dict)

    def add_hit(self, hit: BlastHit) -> None:
        self.hits_by_chrom.setdefault(hit.sseqid, []).append(hit)

    @property
    def chroms(self) -> list[str]:
        return list(self.hits_by_chrom.keys())

    def hits_for(self, chrom: str) -> list[BlastHit]:
        return self.hits_by_chrom.get(chrom, [])

    @property
    def num_chroms(self) -> int:
        return len(self.hits_by_chrom)

    def select_top_n(self, n: int) -> "HitGroup":
        """
        Return a new HitGroup containing only the first N chromosomes
        in the order they appeared in the BLAST file.
        Matches C++ generateIDsFileFromBlastTab selection logic exactly.
        """
        selected = dict(list(self.hits_by_chrom.items())[:n])
        g = HitGroup(query_id=self.query_id)
        g.hits_by_chrom = selected
        return g


def parse_blast_tab(path: str | Path) -> list[HitGroup]:
    """
    Parse a BLAST tab file and return one HitGroup per unique query ID.
    Skips blank lines and comment lines starting with '#'.
    Raises ValueError on malformed lines.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"BLAST file not found: {path}")

    groups: dict[str, HitGroup] = {}
    order:  list[str] = []

    with path.open() as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 12:
                raise ValueError(
                    f"{path}:{lineno}: expected 12 fields, "
                    f"got {len(parts)}"
                )

            try:
                hit = BlastHit(
                    qseqid   = parts[0],
                    sseqid   = parts[1],
                    pident   = float(parts[2]),
                    length   = int(parts[3]),
                    mismatch = int(parts[4]),
                    gapopen  = int(parts[5]),
                    qstart   = int(parts[6]),
                    qend     = int(parts[7]),
                    sstart   = int(parts[8]),
                    send     = int(parts[9]),
                    evalue   = float(parts[10]),
                    bitscore = float(parts[11]),
                )
            except (ValueError, IndexError) as exc:
                raise ValueError(
                    f"{path}:{lineno}: parse error: {exc}"
                ) from exc

            if hit.qseqid not in groups:
                groups[hit.qseqid] = HitGroup(query_id=hit.qseqid)
                order.append(hit.qseqid)
            groups[hit.qseqid].add_hit(hit)

    return [groups[qid] for qid in order]


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python blast.py <blast.tab>")
        sys.exit(1)

    groups = parse_blast_tab(sys.argv[1])
    print(f"Parsed {len(groups)} query group(s):")
    for g in groups:
        total = sum(len(v) for v in g.hits_by_chrom.values())
        print(f"  [{g.query_id}]  "
              f"{g.num_chroms} chromosomes, {total} total hits")
        for chrom, hits in g.hits_by_chrom.items():
            strands = set(h.strand for h in hits)
            coords  = [(h.chrom_start, h.chrom_end) for h in hits]
            print(f"    {chrom}: {len(hits)} hit(s)  "
                  f"strand={strands}  coords={coords}")
