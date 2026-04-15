from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO


# ----------------------------
# Configuration
# ----------------------------

EMAIL = "your_email@example.com"   # replace with your email
TOOL = "biopython_blast_assignment"
DATABASE = "nr"

HBA_FILE = "hba_p69905.fasta"      # protein FASTA
HBB_FILE = "hbb_nm_000518_5.fasta" # nucleotide FASTA transcript

TOP_N = 10
EXCLUDE_HUMAN = True


# ----------------------------
# Data structures
# ----------------------------

@dataclass
class BlastHit:
    target_sequence: str
    species: str
    length: int
    evalue: float
    percent_identity: float


# ----------------------------
# Helpers
# ----------------------------

def read_single_fasta_sequence(path: str) -> str:
    record = SeqIO.read(path, "fasta")
    return str(record.seq)


def species_from_hit_def(hit_def: str) -> str:
    """
    Tries to extract species from BLAST hit_def by taking the last [...] group.
    Example:
      'hemoglobin subunit alpha [Gorilla gorilla gorilla]'
      -> 'Gorilla gorilla gorilla'
    """
    if "[" in hit_def and "]" in hit_def:
        start = hit_def.rfind("[")
        end = hit_def.rfind("]")
        if start != -1 and end != -1 and end > start:
            return hit_def[start + 1:end].strip()
    return "Unknown"


def target_name_from_hit_def(hit_def: str) -> str:
    """
    Removes the trailing [species] from the hit description.
    """
    if "[" in hit_def:
        return hit_def[:hit_def.rfind("[")].strip()
    return hit_def.strip()


def is_human_species(species: str) -> bool:
    return species.lower() == "homo sapiens"


def percent_identity(hsp) -> float:
    if hsp.align_length == 0:
        return 0.0
    return (hsp.identities / hsp.align_length) * 100.0


def pretty_print_table(title: str, hits: list[BlastHit]) -> None:
    print(f"\n{title}")
    print("-" * len(title))

    headers = ["target_sequence", "species", "length", "E-value", "%identity"]
    rows = [
        [
            h.target_sequence,
            h.species,
            str(h.length),
            f"{h.evalue:.3e}",
            f"{h.percent_identity:.2f}",
        ]
        for h in hits
    ]

    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def fmt(row: Iterable[str]) -> str:
        return " | ".join(cell.ljust(widths[i]) for i, cell in enumerate(row))

    print(fmt(headers))
    print("-+-".join("-" * w for w in widths))
    for row in rows:
        print(fmt(row))


# ----------------------------
# BLAST logic
# ----------------------------

def run_remote_blast(program: str, sequence: str, database: str = DATABASE):
    """
    Run BLAST remotely via NCBI.
    program: blastp or blastx
    sequence: FASTA sequence string (no header needed)
    """
    NCBIWWW.email = EMAIL
    NCBIWWW.tool = TOOL

    # XML is the safest format for parsing automatically
    handle = NCBIWWW.qblast(
        program=program,
        database=database,
        sequence=sequence,
        format_type="XML"
    )
    return handle


def parse_top_unique_species(xml_handle, top_n: int = 10) -> list[BlastHit]:
    """
    Parse BLAST XML and keep the first top_n hits from distinct species,
    excluding human if requested.
    """
    record = NCBIXML.read(xml_handle)

    unique_species: set[str] = set()
    selected: list[BlastHit] = []

    for alignment in record.alignments:
        if not alignment.hsps:
            continue

        hsp = alignment.hsps[0]
        species = species_from_hit_def(alignment.hit_def)
        target = target_name_from_hit_def(alignment.hit_def)

        if EXCLUDE_HUMAN and is_human_species(species):
            continue

        if species in unique_species:
            continue

        unique_species.add(species)
        selected.append(
            BlastHit(
                target_sequence=target,
                species=species,
                length=alignment.length,
                evalue=hsp.expect,
                percent_identity=percent_identity(hsp),
            )
        )

        if len(selected) == top_n:
            break

    return selected


def process_query(label: str, fasta_file: str, program: str) -> list[BlastHit]:
    seq = read_single_fasta_sequence(fasta_file)
    xml_handle = run_remote_blast(program=program, sequence=seq)
    hits = parse_top_unique_species(xml_handle, top_n=TOP_N)
    pretty_print_table(label, hits)
    return hits


# ----------------------------
# Main
# ----------------------------

def main() -> None:
    print("Running HBA with BLASTP against nr...")
    hba_hits = process_query(
        label="HBA (P69905) - top non-human protein hits",
        fasta_file=HBA_FILE,
        program="blastp",
    )

    # Be polite to NCBI before a second remote query
    time.sleep(10)

    print("\nRunning HBB transcript with BLASTX against nr...")
    hbb_hits = process_query(
        label="HBB transcript (NM_000518.5) - top non-human protein hits",
        fasta_file=HBB_FILE,
        program="blastx",
    )

    # Example of accessing the top non-human hit E-value
    if hba_hits:
        print(f"\nTop non-human HBA E-value: {hba_hits[0].evalue:.3e}")
    if hbb_hits:
        print(f"Top non-human HBB E-value: {hbb_hits[0].evalue:.3e}")


if __name__ == "__main__":
    main()