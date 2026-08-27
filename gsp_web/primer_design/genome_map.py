"""
genome_map.py
Maps genome display names (as used in BLAST_DATABASES) to their
GFF3 annotation file and chromosome sequence-ID prefix.
"""

import os
from django.conf import settings

# ---------------------------------------------------------------------------
# Each entry: display_name -> {
#   'gff'    : filename inside GFF_DIR (may be HC, LC, or single file),
#   'gff_hc' : HC-only file (optional),
#   'gff_lc' : LC-only file (optional),
#   'chr_prefix': seqid prefix used in that GFF (e.g. 'Ta_Jagger_')
# }
# ---------------------------------------------------------------------------
_MAP = {
    # ── Bread wheat ──────────────────────────────────────────────────────
    # gnl:            the gnl library name stored in the BLAST DB
    # blast_chr_prefix: prefix used for chromosome seq IDs inside the BLAST DB
    #                   (may differ from chr_prefix in the GFF)
    'Chinese Spring (IWGSC ref1 2018)': {
        'gff_hc':          'T_aestivum_CS_ref1_HC.sorted.gff3.gz',
        'gff_lc':          'T_aestivum_CS_ref1_LC.sorted.gff3.gz',
        'chr_prefix':       'Ta_CS1_',
        'gnl':              'IWGSC_rf1',
        'blast_chr_prefix': 'chr',          # BLAST uses chr1A, chr2B, etc.
    },
    'Chinese Spring (IWGSC ref2 2021)': {
        'gff_hc':          'T_aestivum_CS_ref2_HC.sorted.gff3.gz',
        'gff_lc':          'T_aestivum_CS_ref2_LC.sorted.gff3.gz',
        'chr_prefix':       'Ta_CS2_',
        'gnl':              'IWGSC_2.1',
        'blast_chr_prefix': 'Ta_CS_',       # BLAST uses Ta_CS_1A, Ta_CS_2B, etc.
    },
    'Jagger (2019)': {
        'gff':             'T_aestivum_Jagger.sorted.gff3.gz',
        'chr_prefix':       'Ta_Jagger_',
        'gnl':              '10Jagger',
        'blast_chr_prefix': 'Ta_Jagger_',
    },
    'Arina (2019)': {
        'gff':             'T_aestivum_ArinaLrFor.sorted.gff3.gz',
        'chr_prefix':       'Ta_Arina_',
        'gnl':              '10Arina',
        'blast_chr_prefix': 'Ta_Arina_',
    },
    'Julius (2019)': {
        'gff':             'T_aestivum_Julius.sorted.gff3.gz',
        'chr_prefix':       'Ta_Julius_',
        'gnl':              '10Julius',
        'blast_chr_prefix': 'Ta_Julius_',
    },
    'Lancer (2019)': {
        'gff':             'T_aestivum_Lancer.sorted.gff3.gz',
        'chr_prefix':       'Ta_Lancer_',
        'gnl':              '10Lancer',
        'blast_chr_prefix': 'Ta_Lancer_',
    },
    'Landmark (2019)': {
        'gff':             'T_aestivum_Landmark.sorted.gff3.gz',
        'chr_prefix':       'Ta_Landmark_',
        'gnl':              '10Landmark',
        'blast_chr_prefix': 'Ta_Landmark_',
    },
    'Mace (2019)': {
        'gff':             'T_aestivum_Mace.sorted.gff3.gz',
        'chr_prefix':       'Ta_Mace_',
        'gnl':              '10Mace',
        'blast_chr_prefix': 'Ta_Mace_',
    },
    'Mattis (2019)': {
        'gff':             'T_aestivum_Mattis.sorted.gff3.gz',
        'chr_prefix':       'Ta_Mattis_',
        'gnl':              '10Mattis',
        'blast_chr_prefix': 'Ta_Mattis_',
    },
    'Norin61 (2019)': {
        'gff':             'T_aestivum_Norin.sorted.gff3.gz',
        'chr_prefix':       'Ta_Norin_',
        'gnl':              '10Norin',
        'blast_chr_prefix': 'Ta_Norin_',
    },
    'Spelta (2019)': {
        'gff':             'T_spelta.sorted.gff3.gz',
        'chr_prefix':       'T_spelta_',
        'gnl':              '10Spelta',
        'blast_chr_prefix': 'Ta_spelta_',   # note lowercase 's' in BLAST
    },
    'Stanley (2019)': {
        'gff':             'T_aestivum_Stanley.sorted.gff3.gz',
        'chr_prefix':       'Ta_Stanley_',
        'gnl':              '10Stanley',
        'blast_chr_prefix': 'Ta_Stanley_',
    },
    'Kariega (2022)': {
        'gff':             'T_aestivum_Kariega.sorted.gff3.gz',
        'chr_prefix':       'Ta_Kariega_',
        'gnl':              'TaKrga',
        'blast_chr_prefix': 'Ta_Krga_',
    },
    'Fielder (2021)': {
        'gff':             'T_aestivum_Fielder.sorted.gff3.gz',
        'chr_prefix':       'Ta_Fielder_',
        'gnl':              'Ta_Fldr',
        'blast_chr_prefix': 'Ta_Fldr_',
    },
    # ── Einkorn wheat ────────────────────────────────────────────────────
    'T. monococcum TA10622 (2022)': {
        'gff_hc':          'T_monococcum_TA10622_HC.sorted.gff3.gz',
        'gff_lc':          'T_monococcum_TA10622_LC.sorted.gff3.gz',
        'chr_prefix':       'Tm_TA10622_',
        'gnl':              'TA10622',
        'blast_chr_prefix': 'Tm_TA10622_',
    },
    'T. monococcum TA299 (2022)': {
        'gff_hc':          'T_monococcum_TA299_HC.sorted.gff3.gz',
        'gff_lc':          'T_monococcum_TA299_LC.sorted.gff3.gz',
        'chr_prefix':       'Tm_TA299_',
        'gnl':              'TA299',
        'blast_chr_prefix': 'Tm_TA299_',
    },
    # ── Durum wheat ──────────────────────────────────────────────────────
    'Svevo durum (2019)': {
        'gff_hc':          'T_dicoccoides_Svevo_HC.sorted.gff3.gz',
        'gff_lc':          'T_dicoccoides_Svevo_LC.sorted.gff3.gz',
        'chr_prefix':       'Ttd_Svevo_',
        'gnl':              'Td_Svvo',
        'blast_chr_prefix': 'whe_Td_AB_Svevo_',
    },
    'Zavitan durum (2019)': {
        'gff':             'T_dicoccoides_Zaviton.sorted.gff3.gz',
        'chr_prefix':       'Ttd_Zavitan2_',
        'gnl':              'WEWv2',
        'blast_chr_prefix': 'Ttd_Zvtn2_',
    },
    # ── Aegilops ─────────────────────────────────────────────────────────
    'Ae. tauschii AL8-78 (2017)': {
        'gff_hc':          'Ae_tauschii_AL8_78_HC.sorted.gff3.gz',
        'gff_lc':          'Ae_tauschii_AL8_78_LC.sorted.gff3.gz',
        'chr_prefix':       'Aet_AL8-78_',
        'gnl':              'AetAL878',
        'blast_chr_prefix': 'Aegilops_',
    },
    'Ae. tauschii 5.0 Luo (2021)': {
        'gff':             'Ae_tauschii_v5.sorted.gff3.gz',
        'chr_prefix':       'Aet_v5_',
        'gnl':              'Aet_v5',
        'blast_chr_prefix': 'Aet_v5_',
    },
    'Ae. tauschii AY17 (2021)': {
        'gff':             'Ae_tauschii_AY17.sorted.gff3.gz',
        'chr_prefix':       'Aet_AY17_',
        'gnl':              'Aet_AY17',
        'blast_chr_prefix': 'Aet_AY17_',
    },
    'Ae. tauschii AY61 (2021)': {
        'gff':             'Ae_tauschii_AY61.sorted.gff3.gz',
        'chr_prefix':       'Aet_AY61_',
        'gnl':              'Aet_AY61',
        'blast_chr_prefix': 'Aet_AY61_',
    },
    'Ae. tauschii XJ02 (2021)': {
        'gff':             'Ae_tauschii_XJ02.sorted.gff3.gz',
        'chr_prefix':       'Aet_XJ02_',
        'gnl':              'Aet_XJ02',
        'blast_chr_prefix': 'Aet_XJ02_',
    },
    'Ae. longissima Sd (2022)': {
        'gff':             'Ae_longissima.sorted.gff3.gz',
        'chr_prefix':       'Ael_',
        'gnl':              'Ael',
        'blast_chr_prefix': 'Ael_',
    },
    'Ae. speltoides Sb (2022)': {
        'gff':             'Ae_speltoides.sorted.gff3.gz',
        'chr_prefix':       'Aes_',
        'gnl':              'Aes',
        'blast_chr_prefix': 'Aes_',
    },
    'Ae. sharonensis Sd (2022)': {
        'gff':             'Ae_sharonensis.sorted.gff3.gz',
        'chr_prefix':       'Aesh_',
        'gnl':              'Aesh',
        'blast_chr_prefix': 'Aesh_',
    },
}


def get_gff_paths(genome_name: str, confidence: str = 'HC') -> list[str]:
    """
    Return list of absolute GFF3 gz paths for a genome.
    confidence: 'HC', 'LC', or 'both'
    """
    entry = _MAP.get(genome_name)
    if not entry:
        return []

    gff_dir = settings.GFF_DIR
    paths = []

    if 'gff_hc' in entry and 'gff_lc' in entry:
        if confidence in ('HC', 'both'):
            paths.append(os.path.join(gff_dir, entry['gff_hc']))
        if confidence in ('LC', 'both'):
            paths.append(os.path.join(gff_dir, entry['gff_lc']))
    elif 'gff' in entry:
        paths.append(os.path.join(gff_dir, entry['gff']))

    return [p for p in paths if os.path.exists(p)]


def get_chr_prefix(genome_name: str) -> str:
    """Return the chromosome seqid prefix for this genome."""
    entry = _MAP.get(genome_name, {})
    return entry.get('chr_prefix', '')


def list_genomes() -> list[str]:
    """All genome display names with GFF support."""
    return list(_MAP.keys())


def gff_chrom_to_blast_entry(chrom: str, genome_name: str) -> str:
    """
    Convert a GFF chromosome name (e.g. 'Ta_Jagger_2D') to the correct
    blastdbcmd -entry argument (e.g. 'gnl|10Jagger|Ta_Jagger_2D').

    Different genomes use different naming conventions in their BLAST databases.
    This function handles all known mismatches between GFF seqids and BLAST entries.
    """
    entry = _MAP.get(genome_name, {})
    gnl              = entry.get('gnl', '')
    gff_prefix       = entry.get('chr_prefix', '')
    blast_prefix     = entry.get('blast_chr_prefix', gff_prefix)

    # Replace the GFF chromosome prefix with the BLAST chromosome prefix
    if gff_prefix and chrom.startswith(gff_prefix):
        blast_chrom = blast_prefix + chrom[len(gff_prefix):]
    else:
        blast_chrom = chrom

    return f"gnl|{gnl}|{blast_chrom}" if gnl else blast_chrom
