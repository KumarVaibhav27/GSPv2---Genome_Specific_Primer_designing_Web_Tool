# GSP — Genome Specific Primer Design Tool
### Version 2.0

GSP is a web-based tool for designing PCR primers that are specific to a single subgenome in polyploid species. It supports 48 genome assemblies spanning the full ploidy range of the Triticeae tribe: 20 barley (*Hordeum vulgare*) accessions, one rye (*Secale cereale* Lo7) assembly, and 27 wheat and *Aegilops* assemblies covering elite cultivars and diploid progenitors.

---

## What's New in v2.0

Version 1 was a command-line C++ tool that required manual sequence preparation and produced plain CSV output. Version 2 is a complete rewrite as a Django web application with three integrated tools, an interactive alignment viewer, and Excel download support.

| Feature | v1 (C++) | v2 (Django) |
|---|---|---|
| Interface | Command line | Web browser |
| Input | Pre-aligned FASTA files | Paste / upload, auto-aligned |
| BLAST integration | Manual | Automated background job with progress bar |
| KASP marker design | Not available | Full SNP detection + KASP primer output |
| Alignment viewer | Not available | Interactive MSA viewer with primer overlay |
| Output | CSV only | CSV + Excel (multi-sheet by chromosome) |
| Supported genomes | 1 | 48 across 6 species groups |

---

## The Three Tools + AI Assistant

A floating AI chatbot (powered by Claude) is available on every page. It can answer questions about primer design and plant genomics, guide users through the gene-to-primer workflow (finding a gene, extracting its sequence, and launching primer design), and explain GSP results. The chatbot is limited to GSP-related and plant genomics topics.

---

## The Three Tools

### Tab 1 — Paste Your Own Sequences
Design 3′-end specific primers when you already have sequences for each homeolog (e.g. the A, B and D copies of a gene).

- Paste one or multiple sequences in FASTA format
- GSP aligns them using MUSCLE, finds positions unique to each sequence's 3′ end, and runs Primer3 to design primers anchored at those positions
- Results are ranked by 3′-end specificity score
- Click any primer pair to visualise its position in the multiple sequence alignment

### Tab 2 — BLAST Against Genome
Design genome-specific primers from a single query sequence, without needing the homeolog sequences in advance.

- Paste a single FASTA query sequence and select a reference genome
- GSP runs BLASTN against the genome in a background thread, extracts all homeologous hits, aligns them with MUSCLE, and identifies subgenome-unique positions
- A live progress bar tracks the job (large genomes typically complete within a few minutes, depending on genome size and server load)
- Results are grouped by chromosome with per-primer location verification
- Download as multi-sheet Excel (one sheet per chromosome group)

### Tab 3 — KASP Marker Design
Design KASP (Kompetitive Allele Specific PCR) markers for SNP genotyping between two parental lines.

- Paste Parent 1 (Wild Type) and Parent 2 (Mutant) sequences
- GSP aligns them, detects all SNPs, and presents them as selectable chips (positions in 1-based coordinates)
- Select one or more SNP positions and click Design KASP Markers
- For each SNP, GSP outputs:
  - **FAM primer** — allele-specific primer for allele 1 (3′ end at SNP)
  - **HEX primer** — allele-specific primer for allele 2 (3′ end at SNP)
  - **Common primer** — shared reverse (or forward) primer
- The alignment viewer shows primer positions with the allele-specific 3′ end highlighted in blue
- Download all markers as Excel

---

## How to Use

### Tab 1 — Paste Your Own Sequences (step by step)

Use this tab when you already have the sequences for each homeolog (e.g. you have the 3A, 3B, and 3D copies of a gene).

1. Paste two or more homeologous sequences in FASTA format into the text box. Each sequence must have a header line starting with `>` (e.g. `>3A`, `>3B`, `>3D`). You can also upload one or more FASTA files.
2. Optionally adjust the amplicon size range (default 150–1500 bp).
3. Click **Design Primers**.
4. Results appear as tabs — one per sequence. Each tab shows up to 20 primer pairs ranked by specificity.
5. Click any primer pair row to highlight its position in the multiple sequence alignment viewer below.
6. Click **Download Excel** to export all results.

> **Note:** Tab 1 only guarantees specificity relative to the sequences you provide. For full genome-wide specificity, use Tab 2.

---

### Tab 2 — BLAST Against Genome (step by step)

Use this tab when you have a single query sequence and want primers that are specific across the entire reference genome, including all homeologs.

1. Paste a single FASTA sequence (one sequence only) into the text box, or upload a FASTA file.
2. Select the reference genome from the dropdown (e.g. Chinese Spring IWGSC ref2.1).
3. Optionally adjust the amplicon size range (default 150–700 bp).
4. Click **Run BLAST + Design Primers**.
5. A progress bar tracks the job — runtime varies with genome size and server load.
6. Results appear grouped by chromosome. Primers marked **3′-specific** have their 3′ end verified genome-wide.
7. Click **Download Excel** to get a multi-sheet workbook (one sheet per chromosome group).

---

### Tab 3 — KASP Marker Design (step by step)

Use this tab to design KASP genotyping assays for a known SNP between two parental lines.

1. Paste the Parent 1 (Wild Type) sequence into the first box and Parent 2 (Mutant) into the second box. Both must cover the same genomic region.
2. Click **Align & Find SNPs**.
3. SNP positions appear as clickable chips above the alignment viewer. Each chip shows the position and the two alleles (e.g. `pos 142: A/T`).
4. Click one or more SNP chips to select them, then click **Design KASP Markers**.
5. For each selected SNP, GSP outputs a FAM primer, a HEX primer, and a common primer. The allele-specific 3′ end is highlighted in blue in the alignment viewer.
6. Click **Download Excel** to export all markers.

---

### Using the AI Chatbot

The floating flask character in the bottom-right corner is the GSP Assistant. Click it to open the chat window.

The chatbot can:
- Answer questions about primer design, PCR, and plant genomics
- Help you find a gene by name or coordinate and extract its sequence
- Automatically launch primer design once a sequence is ready
- Explain your results

**Example workflow via chatbot:**
1. Type: *"I want primers for TaGW2 in the Jagger genome"*
2. The chatbot searches for the gene, shows a clickable gene list, and asks you to pick one.
3. Click a gene in the table → the chatbot extracts the sequence.
4. Type *"yes"* to confirm → primer design launches automatically in Tab 2.

The chatbot is limited to GSP-related and plant genomics topics.

---

### Interpreting Results

| Column | What it means |
|---|---|
| `3′-end L / R` | **Y** = the primer's 3′ end hits exactly one location genome-wide. This is the most important filter — always prefer Y/Y pairs. |
| `Score` | Combined rank based on Tm balance, GC content, amplicon size, and number of unique positions. Higher is better. |
| `Unique sites L / R` | Number of positions in the MSA where this primer differs from all other homeologs. More is better. |
| `Left/Right Tm` | Melting temperature (°C). Aim for both primers within 2°C of each other. |
| `Product size` | Expected amplicon size in bp. |

### Common errors

**"No specific primer pairs found"** — The sequences may be too similar, or no primer could pass the 3′-end specificity filter. Try a different genomic region (intronic sequences typically have more unique positions than exons).

**"No unique positions found"** — The sequences are nearly identical. Check that you have pasted the correct homeolog sequences, not duplicates of the same one.

**"Sequence is too short"** — Each sequence must be at least 100 bp. For best results, use sequences of at least 500 bp so Primer3 has enough flanking sequence to work with.

---

## How It Works — The Pipeline

### Specificity scoring (Tabs 1 & 2)
GSP uses a two-layer specificity model:

1. **3′-end specificity** — the 3′ end of each primer is BLASTed against the genome. A hit is counted if it has ≥95% identity over ≥18 bp aligned at the primer's 3′ terminus. A primer is 3′-specific if it hits exactly one genomic location under these criteria. This is the primary filter because mismatches at the 3′ end of a primer prevent extension.

2. **Full-length specificity** — the full primer sequence is BLASTed to check for partial off-target binding (≥95% identity, ≥85% of primer length covered). Used as a secondary rank.

Primers passing the 3′-end specificity check are ranked first. Within that group, pairs are ranked by combined score (Tm balance, GC content, amplicon size, position in unique region).

### KASP primer logic
For each selected SNP position, GSP tries two orientations:

- **Forward** — allele-specific primers are forward primers with their 3′ end at the SNP; the common primer is a downstream reverse primer
- **Reverse** — allele-specific primers are reverse primers binding the top strand with their 3′ end at the SNP position; the common primer is an upstream forward primer

The forward orientation is tried first. If Primer3 cannot find a suitable primer in that orientation (due to sequence constraints near the SNP), the reverse orientation is used. The `orientation` field in the output indicates which was selected.

---

## Supported Reference Genomes

48 genome assemblies across 6 species groups:

**Barley (*Hordeum vulgare* H — diploid) — 20 accessions**
Morex v3 (2021), Morex v1 (2017), and 18 additional accessions from the barley pan-genome (Jayakodi et al., 2020; Mascher et al., 2021)

**Rye (*Secale cereale* R — diploid) — 1 assembly**
Lo7 (2021; Rabanus-Wallace et al., 2021)

**Bread wheat (*T. aestivum* ABD — hexaploid) — 14 assemblies**
Chinese Spring IWGSC ref1 (2018), Chinese Spring IWGSC ref2.1 (2021), Kariega (2022), Fielder (2021), Jagger, Arina, Julius, Lancer, Landmark, Mace, Mattis, Norin61, Spelta, Stanley (all 10+ Genome Project 2019)

**Durum wheat (*T. turgidum* AB — tetraploid) — 2 assemblies**
Svevo (2019), Zavitan (2019)

**Einkorn wheat (*T. monococcum* A — diploid) — 2 assemblies**
TA299 (2022), TA10622 (2022)

**Aegilops species — 9 assemblies**
*Ae. tauschii* — AL8-78 (2017), 5.0 Luo (2021), AY17, AY61, T093, XJ02 (2021)
*Ae. longissima* Sd (2022), *Ae. speltoides* Sb (2022), *Ae. sharonensis* Sd (2022)

---

## Project Structure

```
GSP/
├── bin/                        # External binaries (fully bundled — no system installs needed)
│   ├── muscle5                 # MUSCLE wrapper script — calls muscle (v3) with v5-compatible flags
│   ├── muscle                  # MUSCLE v3 binary (Linux x86-64)
│   ├── primer3_core            # Primer3 — primer design engine
│   ├── ntthal                  # Primer3 thermodynamics binary — ΔG calculations
│   └── bedtools                # BEDTools (symlink to Homebrew on macOS)
│
├── blast_db/                   # BLAST nucleotide databases (one folder per genome)
│   └── *.nal                   # Alias files grouping multi-volume databases
│
├── gsp_python/                 # Core Python library (used by the Django views)
│   ├── alignment.py            # MUSCLE wrapper, MSA parsing, unique-position detection
│   ├── blast.py                # BLASTN runner, hit parsing, homeolog extraction
│   ├── primers.py              # Primer3 interface, specificity scoring, ranking
│   ├── kasp.py                 # SNP detection, KASP primer design (FAM/HEX/common)
│   ├── output.py               # Excel and CSV export helpers
│   └── gsp.py                  # Top-level pipeline orchestrator (used by Tab 1 & 2)
│
├── gsp_web/                    # Django web application
│   ├── manage.py
│   ├── db.sqlite3
│   ├── gsp_web/                # Django project settings
│   │   ├── settings.py         # Tool paths, BLAST database registry, ploidy config
│   │   └── urls.py
│   ├── media/jobs/             # Per-job output files (auto-created, auto-populated)
│   └── primer_design/          # Django app
│       ├── views.py            # All three tab endpoints + BLAST background job logic
│       ├── chat_views.py       # AI chatbot endpoint (Claude-powered GSP assistant)
│       ├── gff_lookup.py       # GFF annotation parser for gene search and coordinate lookup
│       ├── genome_map.py       # Genome registry — maps genome names to BLAST DB and GFF paths
│       ├── urls.py             # URL routing
│       └── templates/
│           └── primer_design/
│               └── index.html  # Single-page UI (all three tabs + alignment viewer + chatbot)
│
└── venv/                       # Python virtual environment
```

---

## Example Data

The `examples/` directory contains a worked example for Tab 2 (BLAST Against Genome):

| File | Description |
|---|---|
| `examples/example_input.fas` | Example query sequence (~1,000 bp) in FASTA format |
| `examples/example_output.xlsx` | Expected Excel output from GSP 2.0 for the above input |
| `examples/example_output.csv` | Same results in CSV format |

To reproduce: paste `example_input.fas` into Tab 2, select Chinese Spring IWGSC RefSeq v1.0, and click Run. The output should match `example_output.xlsx`.

---

## Setup & Running

### Requirements
- macOS (Apple Silicon or Intel) or Linux
- Python 3.11+
- BLAST+ (`blastn`, `blastdbcmd`) — install via Homebrew: `brew install blast`
- BEDTools — `brew install bedtools`
- MUSCLE — bundled in `bin/muscle` and `bin/muscle5`
- Primer3 — bundled in `bin/primer3_core`
- **ntthal** (part of Primer3) — required for thermodynamic calculations (ΔG values). Bundled in `bin/ntthal` alongside `primer3_core`. No extra installation needed.

### First-time setup

```bash
cd ~/Desktop/GSP

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate

# Install Python dependencies
pip install django openpyxl anthropic

# ntthal is already bundled in bin/ — no extra step needed
# (If you replace the Primer3 binary, copy the matching ntthal into bin/ as well)
```

### Configuration

All tool paths and settings are controlled via environment variables with sensible macOS defaults built in. No manual editing of `settings.py` is required for a standard macOS setup.

| Environment Variable | Default (macOS) | Description |
|---|---|---|
| `GSP_MUSCLE_PATH` | `bin/muscle5` | Path to MUSCLE binary or wrapper |
| `GSP_PRIMER3_PATH` | `bin/primer3_core` | Path to Primer3 binary |
| `GSP_NTTHAL_PATH` | `bin/ntthal` | Path to ntthal binary (thermodynamics; ships with Primer3) |
| `GSP_BLASTN_PATH` | `/opt/homebrew/bin/blastn` | Path to blastn binary |
| `GSP_BLASTDBCMD_PATH` | `/opt/homebrew/bin/blastdbcmd` | Path to blastdbcmd binary |
| `GSP_BLAST_DB_ROOT` | `../blast_db/` | Root directory containing BLAST databases |
| `GSP_GENOMES_ROOT` | `../genomes/` | Root directory containing reference FASTA files |
| `DJANGO_SECRET_KEY` | insecure default | Set to a long random string in production |
| `DJANGO_DEBUG` | `true` | Set to `false` in production |
| `DJANGO_ALLOWED_HOSTS` | `localhost,127.0.0.1` | Comma-separated list of allowed hostnames |
| `ANTHROPIC_API_KEY` | *(required)* | API key for the Claude-powered chatbot assistant |

### Start the server (macOS)

```bash
cd ~/Desktop/GSP
source venv/bin/activate
python gsp_web/manage.py runserver
```

Open your browser at **http://127.0.0.1:8000**

### Start the server (Linux / custom paths)

Pass environment variables to override the defaults:

```bash
cd /path/to/GSP
source venv/bin/activate

GSP_BLAST_DB_ROOT=/path/to/blast_db \
GSP_PRIMER3_PATH=/path/to/primer3_core \
GSP_MUSCLE_PATH=/path/to/muscle5 \
GSP_BLASTN_PATH=/usr/bin/blastn \
GSP_BLASTDBCMD_PATH=/usr/bin/blastdbcmd \
DJANGO_ALLOWED_HOSTS=your-server-ip,localhost \
DJANGO_DEBUG=false \
DJANGO_SECRET_KEY=your-long-random-secret-key \
ANTHROPIC_API_KEY=your-anthropic-api-key \
python gsp_web/manage.py runserver 0.0.0.0:8000
```

---

## Key Design Decisions

**Why genomic sequences over CDS?**
Introns diverge faster between homeologs than coding exons. A primer designed in an intronic region typically has more unique 3′-end positions available, resulting in higher specificity scores.

**Why MUSCLE v5?**
MUSCLE v5 (`-align` / `-output` flags) replaced the older v3 binary. It is faster and handles longer sequences better. Note: v3 and v5 use different command-line syntax — if you replace the binary, ensure it is v5.

**Why background threads for Tab 2?**
BLAST against a large genome (wheat ~14 GB) can take several minutes depending on genome size and server load. Running this synchronously would time out the HTTP connection. Tab 2 spawns a background thread, writes progress to an in-memory dict (`_JOBS`), and the frontend polls every 2 seconds via `/job-progress/<job_id>/`.

**Why Primer3 `SEQUENCE_FORCE_LEFT_END` for KASP?**
KASP primers must have their 3′ end exactly at the SNP position — this is what makes them allele-specific. `SEQUENCE_FORCE_LEFT_END` forces Primer3 to anchor the primer's 3′ base at the specified position, then designs the rest of the primer upstream.

---

## Output Fields

### Tabs 1 & 2 — Primer pairs
| Field | Description |
|---|---|
| `left_seq` / `right_seq` | Forward and reverse primer sequences (5′→3′) |
| `left_tm` / `right_tm` | Melting temperature (°C) |
| `left_gc` / `right_gc` | GC content (%) |
| `amp_size` | Expected amplicon size (bp) |
| `left_3end_specific` | True if forward primer 3′ end hits exactly 1 genome location |
| `right_3end_specific` | True if reverse primer 3′ end hits exactly 1 genome location |
| `score` | Combined specificity + thermodynamic rank score |

### Tab 3 — KASP markers
| Field | Description |
|---|---|
| `snp_pos` | SNP position in Parent 1 sequence (1-based) |
| `allele1` / `allele2` | Parent 1 (WT) and Parent 2 (Mut) bases at the SNP |
| `fwd1_seq` / `fwd2_seq` | FAM and HEX allele-specific primer sequences |
| `rev_seq` | Common primer sequence |
| `orientation` | `fwd` or `rev` — which strand the allele-specific primers bind |
| `product_size` | Expected amplicon size (bp) |
| `fwd1_tm` / `fwd2_tm` / `rev_tm` | Melting temperatures (°C) |

---
