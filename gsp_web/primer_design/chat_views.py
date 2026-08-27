"""
chat_views.py
Floating chatbot endpoint powered by Claude API with native tool use.

POST /chat/
  Body (JSON):
    {
      "messages": [{"role": "user"|"assistant", "content": "..."}],
      "context":  {
        "active_tab": "tab1"|"tab2"|"tab3",
        "sequence":   "...",
        "database":   "Jagger (2019)",
        "results_summary": "..."
      }
    }

  Returns (JSON):
    {
      "reply": "...",
      "action": null | {
        "type": "run_primer_design",
        "gene_id": "...",
        "sequence": "...",
        "database": "..."
      }
    }
"""

from __future__ import annotations

import json
import logging
import os
from typing import Optional

from django.conf import settings
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_POST

from .genome_map import get_gff_paths, get_chr_prefix, list_genomes, gff_chrom_to_blast_entry
from .gff_lookup import (
    lookup_gene_by_id,
    search_genes_by_keyword,
    search_genes_by_region,
    get_cds_intervals,
    get_intron_intervals,
    extract_sequence,
    list_chromosomes,
    _open_tabix,
    _parse_attrs,
    _revcomp,
)

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# System prompt
# ---------------------------------------------------------------------------

_SYSTEM = """You are GSP Assistant, a helpful AI embedded in the Genome-Specific Primer (GSP) design web tool for plant genomics researchers.

## Tone and style
- Write in plain, conversational prose. No markdown headers (##, **bold**, etc.).
- Keep replies short — 2–4 sentences for simple questions. Only go longer when the user genuinely needs a detailed explanation.
- Never bullet-point everything. If you list things, weave them into a sentence or use a simple comma-separated list.
- One question at a time. Never ask multiple things at once.

## What this tool does
Three tabs, each for a different starting point:

Tab 1 (Paste your own sequences): User pastes 2+ homeologous sequences in FASTA. GSP aligns them, finds discriminating positions, and designs 3′-specific primers. Best when the user already has the sequences.

Tab 2 (BLAST against genome): User pastes one sequence. GSP BLASTs it against a reference genome to find its location, then designs genome-specific primers with 3′-end uniqueness. This is what the chatbot uses when designing from a gene ID.

Tab 3 (KASP markers): User pastes two allele sequences for a SNP. GSP designs KASP assay primers (FAM + HEX allele-specific, plus common reverse) for genotyping.

## Gene-to-Primer workflow
When a user wants primers from a gene:

Step 1 — Find the gene. If they give a Gene ID, call lookup_gene. If not, ask which genome, then which chromosome, then a keyword or coordinate range, then call search_genes. When the tool returns "GENE_LIST_READY:chrom:N", the frontend will automatically display the full gene list as a clickable table. Your reply should ONLY say something like "I found N genes on chromosome X — click any one to select it and I'll extract the sequence." Do NOT list the genes yourself; the table is already shown.

Step 2 — Ask whether they want the CDS or intronic region.

Step 3 — Call extract_gene_sequence and tell the user how many bases were extracted.

Step 4 — Ask if they want to run primer design on that sequence. When they confirm, it is handled automatically.

## Annotation availability
Most wheat genomes in this tool (Jagger, Arina, Lancer, etc.) have structural-only GFF annotations — just gene positions, no functional descriptions. Keyword search only works for Chinese Spring (IWGSC ref2 2021), which has rich Note/Name annotations.

When the search tool returns a result starting with "NO_FUNCTIONAL_ANNOTATIONS:", immediately do a coordinate-based search instead using the known regions below. Do NOT ask the user for permission. Do NOT suggest switching to Chinese Spring. Just search the region, then show every gene returned and ask the user to pick one.

Known regions for common wheat genes:
- Ppd-D1 (chr 2D): 33,000,000–36,000,000
- Ppd-B1 (chr 2B): 33,000,000–36,000,000
- Vrn-A1 (chr 5A): 575,000,000–585,000,000
- Rht-D1 (chr 4D): 495,000,000–505,000,000
- TaGW2 (chr 6A): 608,000,000–616,000,000

For genes not in this list, tell the user keyword search won't work and ask them for approximate coordinates.

## About this tool
GSP (Genome-Specific Primer) was developed by Kumar Vaibhav in the Genomics of Small Grain Cereal Crops lab, Department of Plant Science and Landscape Architecture, University of Maryland, College Park (Maryland, 20742). It is a web-based tool for designing PCR primers that are specific to a single subgenome in polyploid species such as hexaploid wheat (Triticum aestivum). If anyone asks who made the tool, who the creator or developer is, or what the tool is called, answer this directly and accurately.

## Scope — what you can and cannot help with
You are a GSP tool assistant. You can help with:
- Questions about this tool — what it does, who made it, how to use it
- Using any of the three GSP tabs (primer design, BLAST, KASP)
- Plant genomics, wheat biology, homeologs, polyploidy
- Primer design concepts, PCR troubleshooting, BLAST interpretation
- Interpreting GSP results

If a user asks about anything completely unrelated to this tool or plant genomics/biology (e.g. cooking, coding, general knowledge, personal advice), politely decline and redirect. Say something like: "I'm here specifically to help with primer design and the GSP tool — I'm not able to help with that, but feel free to ask me anything about GSP or plant genomics."

## Rules
- Never invent gene IDs, coordinates, or sequences.
- When search_genes returns "GENE_LIST_READY:...", the frontend renders the full table automatically. Just tell the user how many genes were found and to click one. Never list genes yourself.
- Never guess which gene the user wants. Always let the user pick from the rendered list.
- Always end your reply with exactly one line:
  SUGGESTIONS: option1 | option2 | option3
  Choose 2–4 short, natural follow-up options the user would actually say next.

## Available genomes
{genomes}

## Current user context
{context}
"""

# ---------------------------------------------------------------------------
# Tool definitions for Anthropic native tool use
# ---------------------------------------------------------------------------

_TOOLS = [
    {
        "name": "lookup_gene",
        "description": "Look up a gene by its ID in the annotation database. Returns chromosome, coordinates, strand, and functional note.",
        "input_schema": {
            "type": "object",
            "properties": {
                "gene_id": {
                    "type": "string",
                    "description": "The gene ID, e.g. TraesJAG1A01G000100"
                },
                "genome": {
                    "type": "string",
                    "description": "The genome display name, e.g. 'Jagger (2019)'. Infer from the gene ID prefix if possible."
                }
            },
            "required": ["gene_id"]
        }
    },
    {
        "name": "search_genes",
        "description": "Search for genes on a chromosome by keyword (function/description) or coordinate range.",
        "input_schema": {
            "type": "object",
            "properties": {
                "genome": {"type": "string", "description": "Genome display name, e.g. 'Jagger (2019)'"},
                "chrom": {"type": "string", "description": "Chromosome name, e.g. '1A', '2B', '3D'"},
                "keyword": {"type": "string", "description": "Keyword to search in gene function/note (optional)"},
                "start": {"type": "integer", "description": "Start coordinate for region search (optional)"},
                "end": {"type": "integer", "description": "End coordinate for region search (optional)"}
            },
            "required": ["genome", "chrom"]
        }
    },
    {
        "name": "list_chromosomes",
        "description": "List all chromosome names available in a genome's annotation.",
        "input_schema": {
            "type": "object",
            "properties": {
                "genome": {"type": "string", "description": "Genome display name"}
            },
            "required": ["genome"]
        }
    },
    {
        "name": "extract_gene_sequence",
        "description": "Extract the DNA sequence for a gene's CDS or intronic regions from the genome database.",
        "input_schema": {
            "type": "object",
            "properties": {
                "gene_id": {"type": "string", "description": "Gene ID"},
                "genome": {"type": "string", "description": "Genome display name"},
                "region_type": {
                    "type": "string",
                    "enum": ["CDS", "INTRON"],
                    "description": "Which region to extract: CDS (coding sequence) or INTRON"
                }
            },
            "required": ["gene_id", "genome", "region_type"]
        }
    }
]

# ---------------------------------------------------------------------------
# Tool execution
# ---------------------------------------------------------------------------

# Stores extracted sequences pending user confirmation for primer design.
# Keyed by session_key → {gene_id: pending_dict} so different users never
# see each other's pending sequences.
_PENDING_SEQUENCES: dict[str, dict] = {}

# Stores gene search results to be rendered by the frontend as a clickable table.
# Keyed by session_key → {chrom: pending_dict}.
_PENDING_GENE_LISTS: dict[str, dict] = {}


def _infer_genome(gene_id: str, context: dict) -> str:
    """Try to infer genome from gene ID prefix or context."""
    gid = gene_id.upper()
    if 'JAG' in gid:
        return 'Jagger (2019)'
    if 'CS1' in gid or 'TRAESCS1A02' in gid or 'TRAESCS' in gid:
        return 'Chinese Spring (IWGSC ref1 2018)'
    if 'CS2' in gid:
        return 'Chinese Spring (IWGSC ref2 2021)'
    if 'TA10622' in gid.lower() or ('TM.' in gid and '10622' in gid):
        return 'T. monococcum TA10622 (2022)'
    if 'TA299' in gid.lower():
        return 'T. monococcum TA299 (2022)'
    # Fall back to context database
    return context.get('database', '')


def _execute_tool(tool_name: str, tool_input: dict, context: dict, session_key: str = '') -> str:
    genome = tool_input.get('genome') or _infer_genome(tool_input.get('gene_id', ''), context)

    if tool_name == 'lookup_gene':
        gene_id = tool_input.get('gene_id', '').strip()
        gff_paths = get_gff_paths(genome, confidence='HC')
        if not gff_paths:
            gff_paths = get_gff_paths(genome, confidence='both')

        result = lookup_gene_by_id(gene_id, gff_paths)
        if not result:
            # try LC too
            gff_all = get_gff_paths(genome, confidence='both')
            result = lookup_gene_by_id(gene_id, gff_all)
        if not result:
            # try all genomes if genome not identified
            for g in list_genomes():
                paths = get_gff_paths(g, confidence='both')
                result = lookup_gene_by_id(gene_id, paths)
                if result:
                    genome = g
                    break

        if not result:
            return f"Gene '{gene_id}' not found in any annotation database."

        note = result['note'] or 'no functional annotation'
        return (
            f"Gene found:\n"
            f"  ID: {result['gene_id']}\n"
            f"  Genome: {genome}\n"
            f"  Chromosome: {result['chrom']}\n"
            f"  Coordinates: {result['start']}–{result['end']}\n"
            f"  Strand: {result['strand']}\n"
            f"  Function: {note}"
        )

    elif tool_name == 'search_genes':
        chrom_raw = tool_input.get('chrom', '')
        keyword   = tool_input.get('keyword', '')
        start     = tool_input.get('start')
        end       = tool_input.get('end')

        gff_paths = get_gff_paths(genome, confidence='HC')
        prefix    = get_chr_prefix(genome)
        chrom     = chrom_raw if chrom_raw.startswith(prefix) else f"{prefix}{chrom_raw}"

        if keyword:
            results = search_genes_by_keyword(keyword, chrom, gff_paths, max_results=500)
        elif start and end:
            results = search_genes_by_region(chrom, int(start), int(end), gff_paths, max_results=500)
        else:
            return "Please provide either a keyword or a coordinate range (start + end) to search."

        if not results and keyword:
            # Check whether the genome GFF has any functional annotations at all
            # by sampling a few genes and seeing if any have Notes
            has_annotations = False
            for gff_path in gff_paths:
                tbx = _open_tabix(gff_path)
                if tbx is None:
                    continue
                try:
                    sample_count = 0
                    for row in tbx.fetch(chrom):
                        cols = row.split('\t')
                        if len(cols) < 9 or cols[2] != 'gene':
                            continue
                        attrs = _parse_attrs(cols[8])
                        if attrs.get('Note') or attrs.get('Name'):
                            has_annotations = True
                            break
                        sample_count += 1
                        if sample_count > 20:
                            break
                except Exception:
                    pass
                finally:
                    tbx.close()
                if has_annotations:
                    break

            if not has_annotations:
                return (
                    f"NO_FUNCTIONAL_ANNOTATIONS:{genome}:{chrom}:{keyword}"
                )
            return f"No genes with '{keyword}' found on {chrom} in {genome}. Try a different keyword or provide a coordinate range."

        if not results:
            return f"No genes found on {chrom} in {genome}. Check the chromosome name or try a coordinate range."

        # Store the gene list so the frontend can render it as a clickable table
        _PENDING_GENE_LISTS.setdefault(session_key, {})[chrom] = {
            'genes':  results,
            'genome': genome,
            'chrom':  chrom,
            'total':  len(results),
        }
        return f"GENE_LIST_READY:{chrom}:{len(results)}"

    elif tool_name == 'list_chromosomes':
        gff_paths = get_gff_paths(genome, confidence='HC')
        chroms = list_chromosomes(gff_paths)
        if not chroms:
            return f"Could not retrieve chromosome list for {genome}."
        return f"Chromosomes in {genome}:\n" + ", ".join(chroms)

    elif tool_name == 'extract_gene_sequence':
        gene_id     = tool_input.get('gene_id', '').strip()
        region_type = tool_input.get('region_type', 'CDS').upper()

        gff_all = get_gff_paths(genome, confidence='both')
        if not gff_all:
            for g in list_genomes():
                paths = get_gff_paths(g, confidence='both')
                result = lookup_gene_by_id(gene_id, paths)
                if result:
                    genome = g
                    gff_all = paths
                    break

        gene_info = lookup_gene_by_id(gene_id, gff_all)
        if not gene_info:
            return f"Gene '{gene_id}' not found — cannot extract sequence."

        chrom  = gene_info['chrom']
        strand = gene_info['strand']
        blast_db = settings.BLAST_DATABASES.get(genome, '')
        blastdbcmd_path = getattr(settings, 'BLASTDBCMD_PATH', 'blastdbcmd')

        if not blast_db:
            return f"No BLAST database configured for {genome}."

        if region_type == 'CDS':
            intervals = get_cds_intervals(gene_id, chrom, gff_all)
            if not intervals:
                intervals = [(gene_info['start'], gene_info['end'], strand)]
        else:
            intervals = get_intron_intervals(gene_id, chrom, gff_all)
            if not intervals:
                return f"No introns found for {gene_id} — it may be a single-exon gene."

        # For minus-strand genes, CDS intervals are sorted low→high by start,
        # but transcription order is high→low. Reverse so we join in mRNA order.
        if strand == '-':
            intervals = list(reversed(intervals))

        blast_entry = gff_chrom_to_blast_entry(chrom, genome)
        log.info("extract_gene_sequence: chrom=%s blast_entry=%s db=%s", chrom, blast_entry, blast_db)

        segments = []
        for s, e, _st in intervals:
            seg = extract_sequence(blast_entry, s, e, '+', blastdbcmd_path, blast_db)
            if seg:
                segments.append(seg)

        if not segments:
            return (
                f"Sequence extraction failed for {gene_id}. "
                f"The BLAST database for {genome} may not be accessible. "
                f"Check the Django log for the blastdbcmd error."
            )

        full_seq = ''.join(segments)
        # Reverse-complement the full concatenated sequence for minus-strand genes
        if strand == '-':
            full_seq = _revcomp(full_seq)
        _PENDING_SEQUENCES.setdefault(session_key, {})[gene_id] = {
            'sequence': full_seq,
            'genome':   genome,
            'gene_id':  gene_id,
            'region':   region_type,
        }
        return (
            f"Successfully extracted {len(full_seq)} bp of {region_type} sequence "
            f"for {gene_id} ({chrom}, strand {strand}).\n"
            f"The sequence is ready. Shall I run genome-specific primer design against {genome}?"
        )

    return f"Unknown tool: {tool_name}"


# ---------------------------------------------------------------------------
# Check for primer design confirmation
# ---------------------------------------------------------------------------

def _check_for_action(messages: list[dict], session_key: str = '') -> Optional[dict]:
    user_pending = _PENDING_SEQUENCES.get(session_key, {})
    if not user_pending:
        return None
    last_user = next(
        (m['content'] for m in reversed(messages)
         if m['role'] == 'user' and isinstance(m['content'], str)),
        ''
    ).lower()
    confirm_words = ('yes', 'yeah', 'sure', 'go ahead', 'proceed', 'run it', 'do it', 'ok', 'okay', 'run')
    if not any(w in last_user for w in confirm_words):
        return None

    gene_id, pending = next(reversed(user_pending.items()))
    _PENDING_SEQUENCES.get(session_key, {}).pop(gene_id, None)
    return {
        'type':     'run_primer_design',
        'sequence': pending['sequence'],
        'database': pending['genome'],
        'gene_id':  gene_id,
        'amp_min':  150,
        'amp_max':  700,
    }


# ---------------------------------------------------------------------------
# Anthropic client
# ---------------------------------------------------------------------------

_client = None

def _get_client():
    global _client
    if _client is None:
        import anthropic  # type: ignore
        api_key = getattr(settings, 'ANTHROPIC_API_KEY', '') or os.environ.get('ANTHROPIC_API_KEY', '')
        if not api_key:
            raise ValueError("ANTHROPIC_API_KEY is not set")
        _client = anthropic.Anthropic(api_key=api_key)
    return _client


# ---------------------------------------------------------------------------
# Main view — agentic loop with native tool use
# ---------------------------------------------------------------------------

@csrf_exempt
@require_POST
def chat(request):
    try:
        body     = json.loads(request.body)
        messages = body.get('messages', [])
        context  = body.get('context', {})
    except Exception:
        return JsonResponse({'error': 'Invalid JSON'}, status=400)

    if not messages:
        return JsonResponse({'error': 'No messages provided'}, status=400)

    # Ensure session exists so session_key is always available
    if not request.session.session_key:
        request.session.create()
    session_key = request.session.session_key

    # Check if user is confirming a pending primer design (scoped to this session)
    action = _check_for_action(messages, session_key)

    # Build system prompt
    genomes_list = '\n'.join(f'  - {g}' for g in list_genomes())
    context_str  = json.dumps(context, indent=2) if context else 'No context provided.'
    system = _SYSTEM.format(genomes=genomes_list, context=context_str)

    try:
        client = _get_client()
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)

    # Agentic loop: keep calling Claude until it produces a text reply
    # (it may call tools multiple times before giving a final answer)
    loop_messages = list(messages)
    final_reply = ''
    gene_list_action = None
    max_iterations = 5

    for _ in range(max_iterations):
        try:
            response = client.messages.create(
                model='claude-sonnet-4-6',
                max_tokens=1024,
                system=system,
                tools=_TOOLS,
                messages=loop_messages,
            )
        except Exception as e:
            log.error("Anthropic API error: %s", e)
            err_str = str(e).lower()
            if any(k in err_str for k in ('rate_limit', 'quota', 'billing', 'overloaded', '429', '529')):
                return JsonResponse({'error': 'USAGE_LIMIT'}, status=503)
            return JsonResponse({'error': f'AI service error: {e}'}, status=500)

        # Collect text from this response
        text_parts = [b.text for b in response.content if hasattr(b, 'text')]
        tool_uses  = [b for b in response.content if b.type == 'tool_use']

        if not tool_uses:
            # No tool calls — this is the final reply
            final_reply = '\n'.join(text_parts).strip()
            break

        # Execute each tool and build tool_result messages
        # First append the assistant turn with tool_use blocks
        loop_messages.append({'role': 'assistant', 'content': response.content})

        tool_results = []
        for tu in tool_uses:
            result_text = _execute_tool(tu.name, tu.input, context, session_key)
            log.info("Tool %s(%s) → %s", tu.name, tu.input, result_text[:120])

            # Detect gene list ready signal and build structured action
            if result_text.startswith('GENE_LIST_READY:'):
                parts = result_text.split(':')
                chrom = parts[1] if len(parts) > 1 else ''
                pending = _PENDING_GENE_LISTS.get(session_key, {}).pop(chrom, None)
                if pending:
                    gene_list_action = {
                        'type':   'show_gene_list',
                        'genes':  pending['genes'],
                        'genome': pending['genome'],
                        'chrom':  pending['chrom'],
                        'total':  pending['total'],
                    }

            tool_results.append({
                'type':        'tool_result',
                'tool_use_id': tu.id,
                'content':     result_text,
            })

        loop_messages.append({'role': 'user', 'content': tool_results})

    if not final_reply:
        final_reply = 'Sorry, I could not generate a response. Please try again.'

    # Extract SUGGESTIONS line from reply
    suggestions = []
    import re as _re
    sugg_match = _re.search(r'SUGGESTIONS:\s*(.+?)$', final_reply, _re.MULTILINE | _re.IGNORECASE)
    if sugg_match:
        suggestions = [s.strip() for s in sugg_match.group(1).split('|') if s.strip()]
        final_reply = final_reply[:sugg_match.start()].rstrip()

    # Gene list action takes precedence over primer-design confirmation
    if gene_list_action:
        action = gene_list_action

    return JsonResponse({'reply': final_reply, 'action': action, 'suggestions': suggestions})
