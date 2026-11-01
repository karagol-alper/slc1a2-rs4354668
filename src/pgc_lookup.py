#!/usr/bin/env python3
"""
pgc_lookup.py -- rs4354668 lookup + equivalence test in PGC summary statistics.

USAGE
  # 1. resolve the variant (needs network)
  python3 pgc_lookup.py resolve --rsid rs4354668

  # 2. extract from a downloaded sumstats file
  python3 pgc_lookup.py extract --sumstats daner_mdd2025_eur.gz \
      --rsid rs4354668 --chr 11 --pos 35400000 --label "PGC MDD2025 EUR"

  # 3. equivalence test on the extracted row
  python3 pgc_lookup.py tost --beta 0.0012 --se 0.0009 --bound 1.02

  # 4. a priori: how tight a bound can this dataset support?
  python3 pgc_lookup.py precision --ncase 688808 --ncontrol 4364225 --maf 0.41

DEPENDENCIES: numpy, scipy, pandas (+ requests for `resolve`)
"""
import argparse, gzip, sys, json
from math import log, exp, sqrt
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Expected control allele frequency for rs4354668, used as a harmonisation
# sanity check. Source: Rodek 2022 Polish controls (n=462) and gnomAD NFE.
EXPECTED_EAF = {"C": 0.41, "G": 0.41, "A": 0.59, "T": 0.59}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


# ---------------------------------------------------------------------------
def resolve(args):
    """Resolve rsID -> coordinates in both builds; check for merged rsIDs."""
    try:
        import requests
    except ImportError:
        sys.exit("pip install requests")

    print(f"Resolving {args.rsid} ...\n")
    for build, server in [("GRCh38", "https://rest.ensembl.org"),
                          ("GRCh37", "https://grch37.rest.ensembl.org")]:
        url = f"{server}/variation/human/{args.rsid}?content-type=application/json"
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            d = r.json()
        except Exception as e:
            print(f"  {build}: FAILED ({e})")
            continue
        print(f"  {build}:")
        for m in d.get("mappings", []):
            print(f"    chr{m['seq_region_name']}:{m['start']}  "
                  f"alleles {m.get('allele_string')}  strand {m.get('strand')}")
        if d.get("synonyms"):
            print(f"    !! synonyms / merged rsIDs: {d['synonyms']}")
            print("       Search the sumstats for ALL of these, not just the query.")
        print(f"    MAF (1000G): {d.get('MAF')}")
    print("\n  ALWAYS extract by POSITION as well as rsID. Summary-statistics")
    print("  files are inconsistent about which dbSNP build's rsIDs they carry.")


# ---------------------------------------------------------------------------
def sniff_columns(header):
    """Map a sumstats header to canonical names. Handles daner and GWAS-SSF."""
    h = {c.lower(): c for c in header}
    def pick(*cands):
        for c in cands:
            if c in h:
                return h[c]
        # daner frequency columns are FRQ_A_<n> / FRQ_U_<n>
        for k in h:
            if any(k.startswith(c) for c in cands):
                return h[k]
        return None
    return {
        "snp":   pick("snp", "rsid", "variant_id", "markername", "rs_id"),
        "chr":   pick("chr", "chromosome", "#chrom", "chrom"),
        "pos":   pick("bp", "pos", "base_pair_location", "position"),
        "ea":    pick("a1", "effect_allele", "allele1"),
        "oa":    pick("a2", "other_allele", "allele2", "non_effect_allele"),
        "beta":  pick("beta", "effect", "log_odds"),
        "or":    pick("or", "odds_ratio", "oddsratio"),
        "se":    pick("se", "standard_error", "stderr"),
        "p":     pick("p", "p_value", "pval", "pvalue"),
        "eaf":   pick("frq_a_", "eaf", "effect_allele_frequency", "freq", "maf"),
        "frq_u": pick("frq_u_"),
        "info":  pick("info", "imputation_info", "infoscore"),
        "n":     pick("n", "neff", "n_total", "nca"),
    }


def extract(args):
    """Pull the target variant from a sumstats file and harmonise alleles."""
    opener = gzip.open if args.sumstats.endswith((".gz", ".bgz")) else open
    hits = []
    with opener(args.sumstats, "rt") as fh:
        header = fh.readline().rstrip("\n").split()
        if len(header) == 1:
            header = header[0].split("\t")
        cols = sniff_columns(header)
        idx = {k: (header.index(v) if v in header else None)
               for k, v in cols.items() if v}
        print(f"Detected columns: "
              f"{ {k: v for k, v in cols.items() if v} }\n")
        for line in fh:
            f = line.rstrip("\n").split()
            if len(f) == 1:
                f = f[0].split("\t")
            match_rs = (idx.get("snp") is not None
                        and f[idx["snp"]] in (args.rsid, *(args.also or [])))
            match_pos = (args.chr and args.pos and idx.get("chr") is not None
                         and idx.get("pos") is not None
                         and f[idx["chr"]].replace("chr", "") == str(args.chr)
                         and abs(int(f[idx["pos"]]) - args.pos) <= args.window)
            if match_rs or match_pos:
                hits.append(dict(zip(header, f)))

    if not hits:
        print("VARIANT NOT FOUND.")
        print("  Before concluding it is absent, check in order:")
        print("   1. merged/synonym rsIDs   (`resolve` reports these)")
        print("   2. genome build mismatch  (GRCh37 vs GRCh38 position)")
        print("   3. chromosome naming      ('11' vs 'chr11')")
        print("   4. MAF / INFO filtering   (many releases prune INFO<0.6)")
        print("   5. if genuinely absent -> use an LD proxy r2>0.8 via LDlink")
        return

    print(f"{len(hits)} row(s) found in {args.label or args.sumstats}\n")
    for row in hits:
        c = sniff_columns(list(row.keys()))
        ea = row[c["ea"]].upper()
        oa = row[c["oa"]].upper()
        if c["or"] and c["or"] in row:
            beta = log(float(row[c["or"]]))
        else:
            beta = float(row[c["beta"]])
        se = float(row[c["se"]])
        eaf = float(row[c["eaf"]]) if c["eaf"] and c["eaf"] in row else None

        print(f"  raw            : {ea}/{oa}  beta={beta:+.5f}  se={se:.5f}  "
              f"p={row.get(c['p'], 'NA')}")
        if c["info"] and c["info"] in row:
            print(f"  INFO           : {row[c['info']]}")
        if eaf is not None:
            print(f"  effect-allele freq: {eaf:.4f}")

        # -------- harmonisation to the C (low-promoter-activity) allele ------
        note = ""
        if ea == "C":
            b_h = beta
        elif ea == "G":                       # opposite strand, same allele
            b_h, note = beta, " [strand-flipped: G == C]"
        elif ea in ("A", "T"):                # effect allele is the other one
            b_h, note = -beta, " [inverted: effect allele was the A/T allele]"
        else:
            print("  !! unrecognised alleles; harmonise manually"); continue

        print(f"  harmonised to C: beta={b_h:+.5f}  "
              f"OR={exp(b_h):.4f} "
              f"({exp(b_h-1.96*se):.4f}-{exp(b_h+1.96*se):.4f}){note}")

        if eaf is not None:
            exp_f = EXPECTED_EAF.get(ea)
            if exp_f and abs(eaf - exp_f) > 0.15:
                print(f"  !! FREQUENCY MISMATCH: EAF {eaf:.3f} vs expected "
                      f"~{exp_f:.2f} for allele {ea}.")
                print("     Alleles are probably mis-aligned. Do NOT proceed.")
            else:
                print("  frequency check: consistent with expected "
                      f"(~{exp_f:.2f}) -- alleles aligned")
        print()


# ---------------------------------------------------------------------------
def tost(args):
    """Two one-sided tests for equivalence against an OR bound."""
    b, se, d = args.beta, args.se, log(args.bound)
    z_lo = (b + d) / se                      # H0: beta <= -d
    z_hi = (b - d) / se                      # H0: beta >= +d
    p_lo = 1 - stats.norm.cdf(z_lo)
    p_hi = stats.norm.cdf(z_hi)
    p_tost = max(p_lo, p_hi)
    a = args.alpha
    ci90 = (b - stats.norm.ppf(1 - a) * se, b + stats.norm.ppf(1 - a) * se)
    ci95 = (b - 1.96 * se, b + 1.96 * se)

    print("EQUIVALENCE TEST (TOST)")
    print(f"  observed        OR = {exp(b):.4f} "
          f"({exp(ci95[0]):.4f}-{exp(ci95[1]):.4f})")
    print(f"  equivalence bound  = OR {args.bound} "
          f"(|log OR| < {d:.4f})")
    print(f"  {int((1-2*a)*100)}% CI on OR  = "
          f"{exp(ci90[0]):.4f}-{exp(ci90[1]):.4f}")
    print(f"  p(lower one-sided) = {p_lo:.3e}")
    print(f"  p(upper one-sided) = {p_hi:.3e}")
    print(f"  p_TOST             = {p_tost:.3e}")
    equiv = (ci90[0] > -d) and (ci90[1] < d)
    print(f"\n  VERDICT: {'EQUIVALENT' if equiv else 'NOT equivalent'} "
          f"to no effect at the OR {args.bound} bound.")
    if equiv:
        # tightest bound supportable
        tight = exp(max(abs(ci90[0]), abs(ci90[1])))
        print(f"  Tightest bound this dataset supports: OR < {tight:.4f}")
        print("\n  Report as: 'we can exclude an effect larger than "
              f"OR {tight:.3f} (TOST p = {p_tost:.1e})' -- a POSITIVE finding")
        print("  of absence, not merely a failure to reject the null.")


# ---------------------------------------------------------------------------
def precision(args):
    """A priori: what SE, and hence what equivalence bound, is achievable?"""
    n_eff = 4 / (1 / args.ncase + 1 / args.ncontrol)
    se = sqrt(1 / (2 * n_eff * args.maf * (1 - args.maf)))
    print("A PRIORI PRECISION")
    print(f"  cases {args.ncase:,}   controls {args.ncontrol:,}   MAF {args.maf}")
    print(f"  effective N            = {n_eff:,.0f}")
    print(f"  expected SE(log OR)   ~= {se:.5f}")
    for a in (0.05,):
        hw = stats.norm.ppf(1 - a) * se
        print(f"  tightest equivalence bound achievable (alpha={a}): "
              f"OR < {exp(hw):.4f}")
    print(f"\n  For comparison, the published candidate-gene claim is OR 1.57.")
    print(f"  Power to detect OR 1.57 at alpha 5e-8 is essentially 100%.")
    print(f"  If this dataset shows no effect, the claim is refuted, not")
    print(f"  merely unreplicated.")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("resolve"); p.add_argument("--rsid", default="rs4354668")
    p.set_defaults(func=resolve)

    p = sub.add_parser("extract")
    p.add_argument("--sumstats", required=True)
    p.add_argument("--rsid", default="rs4354668")
    p.add_argument("--also", nargs="*", help="synonym rsIDs from `resolve`")
    p.add_argument("--chr", type=int)
    p.add_argument("--pos", type=int)
    p.add_argument("--window", type=int, default=0)
    p.add_argument("--label", default=None)
    p.set_defaults(func=extract)

    p = sub.add_parser("tost")
    p.add_argument("--beta", type=float, required=True)
    p.add_argument("--se", type=float, required=True)
    p.add_argument("--bound", type=float, default=1.02)
    p.add_argument("--alpha", type=float, default=0.05)
    p.set_defaults(func=tost)

    p = sub.add_parser("precision")
    p.add_argument("--ncase", type=int, default=688808)
    p.add_argument("--ncontrol", type=int, default=4364225)
    p.add_argument("--maf", type=float, default=0.41)
    p.set_defaults(func=precision)

    a = ap.parse_args(); a.func(a)
