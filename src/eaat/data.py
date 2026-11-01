"""Data loading for the rs4354668 re-appraisal.

The genotype counts live in exactly one place: data/genotype_counts.csv.
They are transcribed verbatim from the published tables named in that file's
`source_table` column. Nothing here imputes, simulates or adjusts.

STRAND HARMONISATION
--------------------
rs4354668 appears in the literature under three encodings (A/C, T/G and
-181A>C) because studies report opposite strands. The alleles are not
palindromic, so harmonisation is unambiguous:

    A == T   and   C == G

with C/G the minor, low-promoter-activity allele in Europeans. The CSV
stores counts ALREADY harmonised to A/C, and records the published coding
so the transcription can be audited against the source paper.
"""
from __future__ import annotations

import csv
from pathlib import Path

__all__ = ["DATA_DIR", "load_genotypes", "load_ledger", "SHARED_CONTROLS"]

DATA_DIR = Path(__file__).resolve().parents[2] / "data"

#: Studies that report the identical control series. Rodek 2022 and Merk 2019
#: give the same counts (AA 154, AC 234, CC 74; n = 462), the same recruitment
#: site, the same bioethics approval number, the same primers and the same
#: BcnI digest. Treating them as independent would double-count 462 people.
SHARED_CONTROLS = {"rodek22", "merk19"}


def load_genotypes(path: str | Path | None = None) -> dict:
    """Return {study_id: {'case': {...}, 'control': {...}, 'meta': {...}}}.

    Genotype dicts are keyed 'AA', 'AC', 'CC' and harmonised to A/C.
    """
    path = Path(path) if path else DATA_DIR / "genotype_counts.csv"
    out: dict = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            sid = row["study_id"]
            rec = out.setdefault(sid, {"meta": {
                "study": row["study"], "year": int(row["year"]),
                "phenotype": row["phenotype"], "ancestry": row["ancestry"],
                "published_coding": row["published_coding"],
                "source_table": row["source_table"],
                "shared_controls": sid in SHARED_CONTROLS,
            }})
            g = {k: int(row[k]) for k in ("AA", "AC", "CC")}
            if sum(g.values()) != int(row["n"]):
                raise ValueError(
                    f"{sid} {row['group']}: genotypes sum to {sum(g.values())}, "
                    f"but n is {row['n']}")
            rec[row["group"]] = g
    return out


def load_ledger(path: str | Path | None = None) -> list[dict]:
    """Return the replication ledger as a list of dicts, in file order.

    Studies with n_approximate = 1 publish no genotype counts; their sample
    sizes are estimates and they enter only the vote-count and
    excess-significance analyses, never a pooled estimate.
    """
    path = Path(path) if path else DATA_DIR / "study_ledger.csv"
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows.append({
                "study_id": r["study_id"], "study": r["study"],
                "year": int(r["year"]), "phenotype": r["phenotype"],
                "ancestry": r["ancestry"],
                "n_cases": int(r["n_cases"]), "n_controls": int(r["n_controls"]),
                "control_freq_C": float(r["control_freq_C"]),
                "significant": bool(int(r["significant"])),
                "direction": int(r["direction"]),
                "counts_published": bool(int(r["counts_published"])),
                "n_approximate": bool(int(r["n_approximate"])),
            })
    return rows


def label(rec: dict) -> str:
    """Short display label, e.g. 'Rodek et al. 2022 (MDD, EUR_PL)'."""
    m = rec["meta"]
    return f"{m['study']} {m['year']} ({m['phenotype']}, {m['ancestry']})"
