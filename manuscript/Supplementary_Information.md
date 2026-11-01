# Supplementary Information

**Internal replication failure and excess significance in a glial glutamate transporter candidate for depression: a forensic re-appraisal of *SLC1A2* rs4354668**

Alper Karagöl and Taner Karagöl

---

## Contents

| Item | Title |
|---|---|
| Table S1 | Search strategies, yield and exclusions |
| Table S2 | Extracted genotype counts with provenance |
| Table S3 | Per-study estimates under three genetic models, with quality control |
| Table S4 | Statistical power and false-positive report probability |
| Table S5 | Replication ledger and excess-significance calculation |
| Table S6 | Diagnostic discrimination and predictive value |
| Table S7 | Exact tests |
| Table S8 | Bayesian reanalysis |
| Table S9 | Quantitative bias analysis of genotyping error |
| Table S10 | Vibration of effects |
| Table S11 | Simulation of the published record |
| Table S12 | Leave-one-out and cumulative meta-analysis |
| Supplementary Methods S1 | Allele harmonisation across strand encodings |
| Supplementary Methods S2 | Pre-specified analysis plan for the confirmatory study |
| Supplementary Software | Analysis code (five scripts) |

---

## Table S1. Search strategies, yield and exclusions

Eight search strategies were executed. Databases: PubMed/MEDLINE, Google Scholar, and hand-search of the reference lists of all eligible studies and of two reviews of glutamate-transporter genetics in mood disorder.

| # | Strategy | Eligible yield |
|---|---|---|
| 1 | rs4354668 × major depressive disorder, case-control | 2 |
| 2 | Other *SLC1A2* variants × MDD (rs752949, rs1885343, rs4755404) | 1 (already captured by #1) |
| 3 | *SLC1A3* × MDD diagnosis | **0** |
| 4 | *SLC1A3* rs2731880 / rs2269272 × psychiatric case-control | 1 (schizophrenia; not MDD) |
| 5 | Post-stroke depression × *SLC1A2* | **0** |
| 6 | Antidepressant response / treatment resistance × *SLC1A2* | 0 case-control |
| 7 | rs4354668 × schizophrenia, case-control | 4 (contextual) |
| 8 | Suicidal behaviour × *SLC1A2* / *SLC1A3* | 3 (phenotype ineligible) |

### Flow of studies

See Figure 3 for the PRISMA 2020 flow diagram. Title and abstract screening was
conducted iteratively across the eight strategies rather than as a single
indexed pass, so counts are reported from full-text assessment onward; every
assessed record is named below.

| Stage | n |
|---|---|
| Reports assessed for eligibility (full text) | 19 |
| Excluded from primary analysis | 17 |
| — schizophrenia case-control | 5 |
| — bipolar / schizophrenia, other designs | 6 |
| — suicidal behaviour | 2 |
| — quantitative endophenotype | 1 |
| — healthy volunteers, personality trait | 1 |
| — post-mortem expression, no case-control | 1 |
| — Parkinson's disease | 1 |
| **Included in primary meta-analysis** | **2** |

Included: Rodek 2022 (Poland, 161 cases / 462 controls) and Thaweethee 2019
(Thailand, 100 / 100); total 261 cases and 562 controls. Four of the excluded
reports (Merk 2019, Zhang 2015, Nagai 2009, Deng 2004) are retained for the
cross-disorder comparison and the excess-significance analysis. No eligible
study of *SLC1A3* and MDD diagnosis was identified.

### Excluded studies, with reasons

| Study | Reason for exclusion |
|---|---|
| Dallaspezia 2012 | Bipolar disorder; within-patient lithium-response outcome |
| Poletti 2014 (×2) | Bipolar disorder; hippocampal grey-matter endophenotype; no control group |
| Spangaro 2012, 2014 | Schizophrenia; cognitive endophenotype |
| Veldic 2019 | Rapid-cycling bipolar subphenotype; ¹H-MRS glutamate; different variants |
| Fiorentino 2015 | BD/SZ rare-variant burden, not rs4354668 |
| Murphy 2011 | Suicidal behaviour, not MDD diagnosis; variant rs4755404 |
| Choi 2018 | Suicidal behaviour in mixed mood disorders; failed to replicate Murphy |
| Matsumoto 2007 | Healthy volunteers; reward-dependence personality trait |
| Kerdsan 2012 | Post-mortem protein expression; no case-control design |
| Appenzeller 2013 | Parkinson's disease |

---

## Table S2. Extracted genotype counts with provenance

All counts reproduced verbatim from published tables. No imputation, no simulation. Harmonised to the A/C encoding (see Supplementary Methods S1).

| Study | Source table | Group | AA | AC | CC | n |
|---|---|---|---|---|---|---|
| Rodek 2022, *J Clin Med* 11:5914 | Table 2 | MDD | 41 | 71 | 49 | 161 |
| | | Control | 154 | 234 | 74 | 462 |
| Merk 2019, *Neuropsychiatr Dis Treat* 15:989 | Table 1 | Schizophrenia | 129 | 188 | 76 | 393 |
| | | Control | 154 | 234 | 74 | 462 |
| Thaweethee 2019, *Asian Biomed* 12:131 | Table 2 | MDD | 8 | 39 | 53 | 100 |
| | | Control | 7 | 47 | 46 | 100 |

**Note on the shared control series.** The Rodek 2022 and Merk 2019 control rows are identical (AA 154, AC 234, CC 74; n = 462). Both studies report the same recruitment site (Regional Blood Donation Centre, Katowice), the same bioethics approval number (KNW/0022/KB1/34/14), the same published primer pair and the same *Bcn*I restriction digest. These are one control series. Any synthesis treating them as independent would double-count 462 participants.

### Derived allele counts

| Study | Group | A alleles | C alleles | C frequency |
|---|---|---|---|---|
| Rodek 2022 | MDD | 153 | 169 | 0.525 |
| | Control | 542 | 382 | 0.413 |
| Merk 2019 | Schizophrenia | 446 | 340 | 0.433 |
| Thaweethee 2019 | MDD | 55 | 145 | 0.725 |
| | Control | 61 | 139 | 0.695 |

Two-population Wright's *F*ST between the Polish and Thai control series = 0.080.

---

## Table S3. Per-study estimates under three genetic models, with quality control

Risk allele = C, the low-promoter-activity allele.

### Association estimates

| Study | Allelic (C vs A) | Recessive (CC vs AC+AA) | Dominant (CC+AC vs AA) |
|---|---|---|---|
| Rodek 2022 (MDD, Poland) | 1.57 (1.21–2.02) | 2.29 (1.51–3.48) | 1.46 (0.98–2.19) |
| Thaweethee 2019 (MDD, Thailand) | 1.16 (0.75–1.78) | 1.32 (0.76–2.31) | 0.87 (0.30–2.49) |
| Merk 2019 (SZ, Poland) | 1.08 (0.89–1.31) | 1.26 (0.88–1.79) | 1.02 (0.75–1.39) |

### Pooled estimates, MDD only (k = 2)

| Model | Fixed (IV) | Random (DerSimonian–Laird) | Random (Hartung–Knapp) | Q | I² | τ² |
|---|---|---|---|---|---|---|
| Allelic | 1.45 (1.16–1.81) | 1.42 (1.07–1.88) | **1.42 (0.23–8.67)** | 1.40 | 29% | 0.019 |
| Recessive | 1.88 (1.35–2.63) | 1.80 (1.05–3.07) | **1.80 (0.06–57.8)** | 2.40 | 58% | 0.111 |
| Dominant | 1.37 (0.94–2.00) | 1.37 (0.94–2.00) | **1.37 (0.15–12.7)** | 0.83 | 0% | 0.000 |

The Hartung–Knapp interval is the primary random-effects result. At k = 2 the DerSimonian–Laird estimator infers τ² from a single degree of freedom and is severely anticonservative.

### Quality control: Hardy–Weinberg equilibrium and heterozygote deficit

Wright's *F* quantifies the shortfall of heterozygotes relative to Hardy–Weinberg expectation. Bootstrap intervals from 4,000 multinomial draws.

| Series | n | HWE χ² | HWE *p* | *F* | 95% CI |
|---|---|---|---|---|---|
| **Rodek MDD cases (PL)** | 161 | 2.16 | 0.14 | **+0.116** | −0.039 to +0.273 |
| Merk SZ cases (PL) | 393 | 0.26 | 0.61 | +0.026 | −0.075 to +0.121 |
| Rodek/Merk controls (PL) | 462 | 0.91 | 0.34 | −0.044 | −0.136 to +0.048 |
| Thaweethee MDD cases (TH) | 100 | 0.05 | 0.83 | +0.022 | −0.178 to +0.207 |
| Thaweethee controls (TH) | 100 | 1.18 | 0.28 | −0.109 | −0.290 to +0.080 |

Every series passes the conventional Hardy–Weinberg screen. The series producing the single positive association nonetheless carries the largest heterozygote deficit, exceeding the schizophrenia cases genotyped on the same bench with the same enzyme. The interval includes zero; this observation is hypothesis-generating, not probative.

### Case-case comparison

Removing the control series from the analysis entirely:

CC genotype, MDD (Rodek) versus schizophrenia (Merk): **OR 1.82 (1.20–2.77), *p* = 0.005.**

### Sex-by-genotype interaction (recessive model, Rodek 2022 Table 5)

| Stratum | CC / non-CC cases | CC / non-CC controls | OR |
|---|---|---|---|
| Women | 35 / 81 | 34 / 204 | 2.59 (1.51–4.44) |
| Men | 14 / 31 | 40 / 184 | 2.08 (1.01–4.26) |

Ratio of odds ratios 1.25 (0.51–3.06); interaction z = 0.48, **p = 0.63**. No evidence of sex modification.

---

## Table S4. Statistical power and false-positive report probability

### Power of the published studies (allelic model, α = 0.05)

| True OR | Rodek (161/462) | Thaweethee (100/100) | Merk (393/462) |
|---|---|---|---|
| 1.02 | 3.5% | 3.1% | 3.9% |
| 1.05 | 5.7% | 4.1% | 7.2% |
| 1.10 | 11.0% | 6.4% | 16.1% |
| 1.20 | 28.9% | 12.8% | 46.1% |
| 1.57 | 93.5% | 50.4% | 99.6% |

### Sample size for 80% power, MAF 0.41, 1:1 design

| True OR | α = 0.05 | α = 5 × 10⁻⁸ |
|---|---|---|
| 1.05 | 13,574 cases | 68,488 cases |
| 1.10 | 3,549 cases | 17,884 cases |
| 1.20 | 967 cases | 4,860 cases |

### False-positive report probability (Wacholder et al. 2004)

Observed *p* ≈ 0.001 for the recessive finding. Power evaluated at the plausible true effect, not the observed one.

| Prior probability | True OR 1.05 | True OR 1.20 | True OR 2.38 |
|---|---|---|---|
| 0.25 | 0.626 | 0.091 | 0.003 |
| 0.05 | 0.914 | 0.388 | 0.019 |
| 0.01 | 0.982 | 0.768 | 0.090 |
| 0.001 | 0.998 | 0.971 | 0.500 |

---

## Table S5. Replication ledger and excess-significance calculation

### Complete ledger of located case-control tests of rs4354668

| Study | Phenotype | Ancestry | Cases | Controls | Significant | Counts published |
|---|---|---|---|---|---|---|
| Rodek 2022 | MDD | European (Poland) | 161 | 462* | Yes | Yes |
| Thaweethee 2019 | MDD | SE Asian (Thailand) | 100 | 100 | No | Yes |
| Merk 2019 | Schizophrenia | European (Poland) | 393 | 462* | No | Yes |
| Zhang 2015 | Schizophrenia | E Asian (China) | ~1,117 | ~2,289 | Yes | No |
| Nagai 2009 | Schizophrenia | E Asian (Japan) | ~200 | ~200 | No | No |
| Deng 2004 | Schizophrenia | E Asian (Japan) | ~100 | ~100 | No | No |

\* Shared control series. Sample sizes for Nagai 2009 and Deng 2004 are approximate and enter only the power-based expected-significance calculation.

Totals: 2,071 cases / 3,613 controls across all psychiatric phenotypes; 261 cases / 562 controls for MDD.

### Sign test

Four studies report a directional estimate; all four point the same way. Exact binomial *p* = 0.125. With k = 4 this test cannot reach significance even at 4/4, and a literature reporting only confirmatory analyses produces directional consistency by construction.

### Excess significance (Ioannidis & Trikalinos 2007)

Observed significant findings O = 2 of 6.

| Assumed true OR | Expected E | O − E | Binomial *p* | Interpretation |
|---|---|---|---|---|
| 1.02 | 0.23 | +1.77 | **0.019** | Excess |
| 1.05 | 0.41 | +1.59 | **0.059** | Excess |
| 1.10 | 0.93 | +1.07 | 0.237 | Consistent |
| 1.20 | 2.19 | −0.19 | 0.707 | Under-reported |
| 1.45 | 4.29 | −2.29 | 0.991 | Under-reported |

### Fragility of the pooled estimate

Number of additional null studies required to drag the pooled MDD fixed-effect estimate below an equivalence bound of OR 1.05:

| Size of each null study | Number required |
|---|---|
| 200 cases | 22 |
| 500 cases | 9 |
| 1,000 cases | 5 |

---

## Table S6. Diagnostic discrimination and predictive value

### Area under the ROC curve, single variant, log-additive model, risk-allele frequency 0.41

| Per-allele OR | AUC |
|---|---|
| 1.02 | 0.504 |
| 1.05 | 0.509 |
| 1.10 | 0.518 |
| 1.20 | 0.534 |
| **1.57 (published estimate)** | **0.584** |
| 2.00 | 0.627 |
| 3.00 | 0.695 |

Chance = 0.500. Screening instruments are generally considered unusable below approximately 0.75.

### Predictive value of the CC genotype

Lifetime prevalence 0.20; CC frequency 0.17.

| Assumed recessive OR | Risk if CC | Risk if not CC | Sensitivity | Specificity |
|---|---|---|---|---|
| 2.38 (as published) | 38.6% | 16.2% | 32.5% | 87.1% |
| 1.50 | 27.7% | 18.4% | 23.3% | 84.8% |
| 1.20 | 23.2% | 19.3% | 19.5% | 83.9% |
| 1.05 | 20.8% | 19.8% | 17.5% | 83.4% |

---

## Table S7. Exact tests

Cell counts are small; asymptotic *p*-values are reported alongside Fisher's exact test and a mid-*p* correction.

| Contrast | χ² *p* | Fisher exact *p* | mid-*p* |
|---|---|---|---|
| Allelic (C vs A) | 5.27 × 10<sup>-4</sup> | 5.56 × 10<sup>-4</sup> | 5.56 × 10<sup>-4</sup> |
| Recessive (CC vs AC+AA) | 7.57 × 10<sup>-5</sup> | 1.35 × 10<sup>-4</sup> | 1.37 × 10<sup>-4</sup> |
| Dominant (CC+AC vs AA) | 6.38 × 10<sup>-2</sup> | 7.55 × 10<sup>-2</sup> | 6.27 × 10<sup>-2</sup> |

The exact tests do not overturn the asymptotic results. The allelic test treats 322 alleles from 161 individuals as independent observations; the genotypic contrasts do not, and are the appropriate basis for inference.

---

## Table S8. Bayesian reanalysis

Observed allelic estimate: log OR = 0.4493, SE = 0.1301 (OR 1.567). Normal prior on the log odds ratio centred at zero.

| Prior SD | Prior 95% OR range | Posterior OR | Shrinkage | BF₁₀ |
|---|---|---|---|---|
| **0.05** | **0.91–1.10** | **1.060** | **87%** | **2.01** |
| 0.10 | 0.82–1.22 | 1.182 | 63% | 7.28 |
| 0.20 | 0.68–1.48 | 1.371 | 30% | 36.12 |
| 0.50 | 0.38–2.66 | 1.523 | 6% | 67.33 |

A prior standard deviation of 0.05 is the one consistent with common-variant effects in depression. Under it the Bayes factor is 2.0, conventionally "barely worth a bare mention," despite *p* < 0.001. The larger priors make the dependence explicit: a Bayes factor above 50 requires believing in advance that effects up to OR 2.7 are plausible for a common variant.

---

## Table S9. Quantitative bias analysis of genotyping error

**Model.** The true case genotype distribution equals the control distribution (null). Each true A allele is mis-scored as C with probability ε, as would occur under over-digestion in PCR-RFLP — the direction required to inflate CC among cases.

**Required rate: ε = 0.245 (24.5%).**

| Genotype | Observed (cases) | Fitted at ε = 0.245 | True distribution under the null |
|---|---|---|---|
| AA | 41 | 30.6 | 53.7 |
| AC | 71 | 81.4 | 81.5 |
| CC | 49 | 49.0 | 25.8 |

Goodness of fit χ² = 4.9 (2 df), *p* = 0.087. The model is fitted to CC and then under-predicts AA by a third: no single mis-scoring rate reproduces the observed distribution.

### Apparent effect generated by plausible error rates

| Per-allele error rate | Apparent recessive OR from a true null |
|---|---|
| 1% | 1.04 |
| 2% | 1.08 |
| 5% | 1.20 |
| 10% | 1.43 |
| — | **observed: 2.29** |

The authors report re-genotyping 5% of samples with 100% concordance. An error rate of 24.5% would have produced approximately 99% discordance across those samples.

**Conclusion.** Differential genotyping error at any plausible rate cannot generate the observed effect. The heterozygote deficit in Table S3 is not a sufficient explanation and is not advanced as one.

---

## Table S10. Vibration of effects

Every odds ratio obtainable from the Rodek 2022 dataset across five inheritance models and three strata.

| Stratum | Contrast | OR | *p* |
|---|---|---|---|
| women | recessive | 2.59 | 0.0005 \* |
| men | CC vs AA | 2.58 | 0.0341 \* |
| women | CC vs AA | 2.50 | 0.0046 \* |
| all | CC vs AA | 2.49 | 0.0003 \* |
| all | recessive | 2.29 | 0.0001 \* |
| men | recessive | 2.08 | 0.0459 \* |
| men | dominant | 1.75 | 0.1340 |
| men | allelic | 1.65 | 0.0300 \* |
| all | allelic | 1.57 | 0.0006 \* |
| women | allelic | 1.52 | 0.0096 \* |
| all | dominant | 1.46 | 0.0646 |
| women | dominant | 1.27 | 0.3502 |
| men | overdominant | 0.94 | 0.8501 |
| all | overdominant | 0.77 | 0.1527 |
| women | overdominant | 0.64 | 0.0510 |

Range 0.64 to 2.59, a 4.0-fold spread; 9 of 15 reach *p* < 0.05 unadjusted. The reported result is not the largest available.

---

## Table S11. Simulation of the published record

Monte Carlo, 40,000 replicates. Publication bias modelled as suppression of non-significant findings; *q* is the probability that a non-significant study is published. Cells give P(≥ 2 of 6 published studies significant). Observed: 2 of 6.

| True OR | *q* = 1.0 | *q* = 0.5 | *q* = 0.2 | *q* = 0.1 |
|---|---|---|---|---|
| 1.00 | 0.9% | 2.9% | 14.2% | 35.4% |
| 1.02 | 1.9% | 6.5% | 25.7% | 53.3% |
| 1.05 | 5.8% | 17.1% | 48.7% | 76.9% |
| 1.10 | 23.4% | 46.9% | 80.0% | 94.5% |
| 1.20 | 70.4% | 87.6% | 98.1% | 99.7% |

The observed record is consistent with a true null plus 80–90% suppression of negative findings, and with a small true effect (OR 1.05–1.10) under complete reporting. It is not consistent with the published effect size.

---

## Table S12. Leave-one-out and cumulative meta-analysis

| Analysis | Result |
|---|---|
| Omitting Rodek et al. 2022 | OR 1.16 (0.75–1.78) |
| Omitting Thaweethee et al. 2019 | OR 1.57 (1.21–2.02) |
| Cumulative to 2019 (k = 1) | OR 1.16 (0.75–1.78) |
| Cumulative to 2022 (k = 2) | OR 1.45 (1.16–1.81) |

At k = 2 these are reported for completeness, not inference: neither can detect the drift a cumulative analysis is designed to reveal. The pooled estimate rests entirely on one of the two studies.

---

## Supplementary Methods S1. Allele harmonisation across strand encodings

rs4354668 appears in the literature under three encodings, reflecting opposite-strand reporting:

| Encoding | Used by | Low-activity allele |
|---|---|---|
| A/C | Rodek 2022; Mallolas 2006 | C |
| T/G | Thaweethee 2019; Merk 2019; Poletti 2014; Zhang 2015 | G |
| −181A>C | functional promoter literature | C |

The alleles are not palindromic, so harmonisation is unambiguous: **A ≡ T and C ≡ G**, with C/G the minor, low-promoter-activity allele in European populations. All counts in this work are harmonised to A/C.

This is not a cosmetic issue. Two published narrative reviews describe the risk allele inconsistently as a consequence of the dual encoding. Any future analysis drawing on these primary studies must align the effect allele before pooling; taking a published effect estimate at face value without alignment inverts its sign, and an inverted null is indistinguishable from a correctly oriented one.

**Verification.** The harmonisation routine in `pgc_lookup.py` was tested on the same underlying result encoded two ways — as `A1 = A, OR = 0.9988` in daner format and as `effect_allele = C, beta = +0.0012` in GWAS-SSF format — and returns an identical harmonised estimate from both. A frequency check against the expected European C-allele frequency (≈0.41) halts the routine if the reported effect-allele frequency diverges by more than 0.15.

---

## Supplementary Methods S2. Pre-specified analysis plan for the confirmatory study

The analyses below are specified in advance and have not been performed. They are reported here so that the plan is on record before the data are examined.

### Tier 1 — variant-level resolution

1. **Brain eQTL lookup.** rs4354668 against *SLC1A2* expression in GTEx v8 brain tissues, PsychENCODE DLPFC, and astrocyte-specific single-cell eQTL panels. This is listed first because if the variant does not regulate the gene in brain, the association question is moot regardless of outcome. The published functional evidence is a luciferase reporter construct in cell culture, not *in vivo* brain expression.
2. **Variant lookup with equivalence testing.** rs4354668 in PGC MDD 2025 (European excluding 23andMe; clinical/interview-ascertained subset; multi-ancestry), PGC3 SCZ and PGC BD summary statistics. Two one-sided tests against an equivalence bound of OR 1.02, reporting the tightest bound supported. A positive finding of absence is the target, not a failure to reject.
3. **Regional extraction.** All variants across *SLC1A2* and *SLC1A3* ± 100 kb, ancestry-stratified, including **rs12294045** — the one *SLC1A2* variant that survived multiple-testing correction in the only study to apply it, and which has not been followed up.

### Tier 2 — gene and gene-set level

MAGMA gene-based association for *SLC1A2* and *SLC1A3* with an LD-matched 1000 Genomes reference, followed by a competitive gene-set test against the pre-specified astrocytic glutamate–glutamine cycle set below, with gene-size, gene-density and LD covariates. Stratified LD-score regression partitioning depression heritability into astrocyte-specific regulatory annotations, conditioned on the full baseline model.

**Pre-specified gene set** (fixed before analysis):

> *SLC1A2*, *SLC1A3*, *SLC1A6*, *SLC1A7*, *GLUL*, *GLS*, *GLS2*, *SLC38A1*, *SLC38A2*, *SLC38A3*, *SLC7A11*, *AQP4*, *GJA1*, *GJB6*, *KCNJ10*, *ALDH1L1*

### Tier 3 — causal inference

TWAS / S-PrediXcan on genetically predicted EAAT1 and EAAT2 expression using GTEx v8 brain models and astrocyte-specific panels. Two-sample Mendelian randomisation with cis-eQTLs as exposure instruments, with MR-Egger and weighted-median sensitivity analyses, Steiger filtering for direction, and colocalisation to distinguish a shared causal variant from linkage-disequilibrium contamination. Instrument strength to be reported as F-statistics; a single-gene cis-MR has limited power and a null is to be treated as uninformative rather than negative. The reverse direction — depression liability to transporter expression — is to be tested in parallel, since the post-mortem literature cannot distinguish cause from consequence.

### Interpretation of a null result

A null across Tiers 1–3 would establish that astrocytic glutamate-clearance *genetics* do not contribute detectably to depression liability, while leaving the post-mortem and pharmacological evidence for astrocytic dysfunction intact. That would implicate acquired, stress-induced or epigenetic dysregulation of EAAT2 rather than inherited variation, and would strengthen rather than weaken the rationale for EAAT2-targeting therapeutics, since a drug target need not be a heritable risk locus.

---

## Supplementary Software

Five scripts, each with the source genotype counts inline and the originating publication named in the header comment. Python 3.13; numpy 2.4; scipy 1.17; matplotlib 3.10. No other dependencies.

| File | Produces |
|---|---|
| `eaat_mdd_meta.py` | Table S2 derived counts; Table S3 association and pooled estimates; Table S4 power, sample size and FPRP; *F*ST |
| `eaat_mdd_meta_v2.py` | Table S2 shared-control detection; Table S3 quality control, case-case test, sex interaction; multiplicity accounting |
| `eaat_ledger_v3.py` | Table S5 ledger, sign test, excess significance, fail-safe N |
| `diagnostic_utility.py` | Table S6 AUC and predictive value |
| `06_extended_analyses.py` | Tables S7–S12 |
| `05_figure1.py` | Figure 1 (PNG at 400 dpi and vector PDF) |
| `pgc_lookup.py` | Supplementary Methods S1 harmonisation routine; Tier 1 lookup and equivalence-testing utility |

Every numerical value reported in the main text and in this Supplementary Information is regenerated by running these scripts; none is transcribed by hand. The accompanying test suite (`tests/test_reproducibility.py`, 44 assertions) fails if any reported value changes.
