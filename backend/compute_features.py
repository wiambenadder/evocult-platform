"""Turn a species' ligand/receptor FASTA into the model's feature row.

Prototype status:
  * X1  (sequence identity)  — computed here, for real, per protein.
  * X2-X5 (structure, docking, biophysical, evolutionary) — STUBBED. These require
    the heavy pipeline (AlphaFold structures, PRODIGY docking, PAML/HyPhy, expression
    data) and are left as NaN so the model imputes them. Wire each TODO to the
    corresponding pipeline stage to move from prototype to full scoring.

The model tolerates missing features (median imputation), so a FASTA-only request
returns a real X1-driven prediction today, and improves monotonically as the other
axes are filled in.
"""
import numpy as np

# Canonical proteins whose X1_<NAME> columns the model consumes.
PROTEIN_COLUMNS = ["LRP6", "FZD8", "EGFR", "LGR5", "BMPR1A", "BMPR1B",
                   "WNT3A", "EGF", "RSPO1", "NOG", "CHRD", "DKK1", "Fc"]

ALIASES = {
    "WNT3A": ["wnt3a", "wnt-3a", "wnt3"], "RSPO1": ["rspo1", "rspondin1", "r-spondin1"],
    "EGF": ["egf", "epidermal growth factor"], "NOG": ["nog", "noggin"],
    "FZD8": ["fzd8", "frizzled8", "frizzled-8"], "LGR5": ["lgr5"],
    "EGFR": ["egfr", "egf receptor", "erbb1"], "BMPR1A": ["bmpr1a", "alk3"],
    "BMPR1B": ["bmpr1b", "alk6"], "LRP6": ["lrp6"], "CHRD": ["chrd", "chordin"],
    "DKK1": ["dkk1", "dickkopf1"], "Fc": ["fc"],
}

# TODO: replace with the curated human UniProt interface reference panel.
HUMAN_REFERENCE = {}  # {"WNT3A": "MAPLGY...", ...}


def parse_fasta(text):
    records, header, seq = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:].strip(), []
        else:
            seq.append("".join(ch for ch in line.upper() if ch.isalpha()))
    if header is not None:
        records.append((header, "".join(seq)))
    return [(h, s) for h, s in records if s]


def match_protein(header):
    h = header.lower()
    for key, al in ALIASES.items():
        if any(a in h for a in al):
            return key
    return None


def identity(query, ref, cap=500):
    """Global-alignment % identity, normalised over the shorter sequence."""
    a, b = query[:cap], ref[:cap]
    n, m = len(a), len(b)
    if not n or not m:
        return np.nan
    GAP, MATCH, MIS = -1, 1, -1
    prev = np.arange(m + 1) * GAP
    prevM = np.zeros(m + 1, dtype=int)
    for i in range(1, n + 1):
        curr = np.empty(m + 1, dtype=int); currM = np.empty(m + 1, dtype=int)
        curr[0] = i * GAP; currM[0] = 0
        for j in range(1, m + 1):
            eq = a[i - 1] == b[j - 1]
            best = prev[j - 1] + (MATCH if eq else MIS); bm = prevM[j - 1] + (1 if eq else 0)
            if prev[j] + GAP > best:
                best = prev[j] + GAP; bm = prevM[j]
            if curr[j - 1] + GAP > best:
                best = curr[j - 1] + GAP; bm = currM[j - 1]
            curr[j] = best; currM[j] = bm
        prev, prevM = curr, currM
    return max(0.0, min(1.0, prevM[m] / min(n, m)))


def features_from_fasta(fasta_text, reference=None):
    """Return {feature_name: value or np.nan} for the model. Only X1_* are filled."""
    reference = reference or HUMAN_REFERENCE
    records = parse_fasta(fasta_text)
    row = {}

    # ---- X1: sequence identity per protein ----
    for header, seq in records:
        key = match_protein(header)
        if key and reference.get(key):
            row[f"X1_{key}"] = identity(seq, reference[key])

    # ---- X2: structural similarity (iRMSD from AlphaFold) ----
    # TODO: fold ortholog + reference, compute interface iRMSD -> X2_* columns.
    # ---- X3: docking energetics (PRODIGY dG, contact maps) ----
    # TODO: dock complexes, compute dG / contacts -> X3_* columns.
    # ---- X4: biophysical / integrative ----
    # TODO: interface hydrophobicity match, delta pI, complex stability -> X4_* columns.
    # ---- X5: evolutionary context ----
    # TODO: dN/dS (PAML/HyPhy), receptor expression -> X5_* columns.

    return row
