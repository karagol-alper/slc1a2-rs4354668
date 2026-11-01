# Internal replication failure and excess significance in a glial glutamate transporter candidate for depression: a forensic re-appraisal of *SLC1A2* rs4354668

**Running title:** Anatomy of a post-GWAS candidate-gene finding

**Authors:** Alper Karagöl^1,2,3,†,* and Taner Karagöl^1,4,†,*

^1 Istanbul University, Istanbul Medical Faculty, Istanbul, Turkey
^2 Present address: Broad Institute of MIT and Harvard, Cambridge, MA, USA
^3 Present address: Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA
^4 Present address: McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, MA, USA

^† These authors contributed equally.
^* Joint corresponding authors.

Alper Karagöl — akaragol@mit.edu — ORCID 0009-0001-7864-0732
Taner Karagöl — karagolt@mit.edu — ORCID 0009-0005-1011-7661

**Word count:** Abstract 248; main text ~4,400
**Display items:** 4 tables, 3 figures; Supplementary Information (Tables S1–S6, Supplementary Software)
**Keywords:** *SLC1A2*; EAAT2; major depressive disorder; candidate gene; excess significance; meta-research; astrocyte

---

## Abstract

**Background.** Impaired astrocytic glutamate clearance is a leading hypothesis in the pathophysiology of major depressive disorder (MDD). The functional *SLC1A2* promoter variant rs4354668 (−181 A>C), which reduces reporter promoter activity by approximately 30%, has been proposed as a susceptibility locus and candidate biomarker. We conducted a systematic re-appraisal of the underlying evidence.

**Methods.** We systematically searched for case-control association studies of *SLC1A2* and *SLC1A3* variants in MDD, extracted genotype counts, and conducted random-effects meta-analysis with small-*k* variance correction. We then applied a pre-specified six-part forensic audit: detection of overlapping control series, genotyping quality control via heterozygote deficit, within-laboratory cross-disorder replication, direct case-case testing, formal interaction testing of reported subgroup effects, and the Ioannidis–Trikalinos test of excess significance.

**Results.** Two eligible MDD studies were identified (261 cases, 562 controls); none existed for *SLC1A3*. Fixed-effect pooling yielded OR 1.45 (95% CI 1.16–1.81), but the Hartung–Knapp interval appropriate at *k* = 2 spanned 0.23–8.67. The single positive study shared its entire 462-participant control series with a previously published schizophrenia study from the same laboratory in which the variant was null in 2.4-fold more patients (OR 1.08, 95% CI 0.89–1.31). The reported female-specific effect showed no interaction (ratio of ORs 1.25, 95% CI 0.51–3.06; *p* = 0.63). Across six studies, two significant findings were observed against 0.23–0.41 expected at plausible effect sizes (excess significance *p* = 0.019–0.059). Under a prior calibrated to known common-variant effects the published estimate shrinks to OR 1.06 with a Bayes factor of 2.0. Differential genotyping error was modelled and found insufficient, requiring an implausible 24.5% mis-scoring rate. Fifteen defensible analyses of the same data span a four-fold range of odds ratios. Single-variant discrimination is capped at AUC 0.584.

**Conclusions.** The *SLC1A2* depression association is unsupported and shows the signature of selective reporting. We provide a six-part audit, executable from published tables alone, for triaging comparable claims. The six-part audit is generalisable to other candidate-gene claims persisting into the post-GWAS era.

---

## Significance

- Candidate-gene associations continue to be published in psychiatry despite well-powered genome-wide evidence against them, and reviews propagate them as established.
- We show that the leading glial glutamate transporter candidate for depression rests on 261 patients worldwide, fails internal replication within the originating laboratory, and produces more significant findings than the literature's own statistical power can support.
- The audit framework we apply — six checks executable from published tables alone, without access to raw data — offers a reusable tool for triaging legacy candidate-gene claims before they are cited further.

---

## 1. Introduction

Astrocytic pathology is among the more robust neuropathological findings in mood disorder. Glial density is reduced in the prefrontal cortex of individuals who died with depression,¹ and expression of *SLC1A2*, encoding the excitatory amino acid transporter EAAT2 that mediates the large majority of cerebral glutamate uptake, is reduced in anterior cingulate and dorsolateral prefrontal cortex,² hippocampus,³ and locus coeruleus astrocytes.⁴ The relationship is not merely correlational: pharmacological blockade of astrocytic glutamate uptake is sufficient to induce anhedonia and anxiety-like behaviour in rodents.⁵–⁶ Modern genome-wide data converge at the cell-type level, with depression heritability enriched in astrocyte and oligodendrocyte lineages.⁷

Against this background, the promoter variant rs4354668 has attracted sustained attention. It is one of a small minority of psychiatric candidate variants with direct functional validation: the C (equivalently G, on the opposite strand) allele abolishes an AP-2 activator binding site, creates a binding site for the repressor GC-binding factor 2, and reduces luciferase reporter promoter activity by approximately 30%.⁸ A 2022 case-control study reported a two-fold excess of the CC genotype in Polish patients with MDD,⁹ and the finding has since entered reviews as evidence that glial glutamate clearance genotype confers depression liability and may have diagnostic value.

This inference sits uneasily alongside two facts. First, well-powered analyses have found no support for any historical depression candidate gene: Border and colleagues tested the 18 most-studied candidates in up to 620,000 individuals and found them no more associated with depression than randomly chosen genes.¹⁰ Second, the current Psychiatric Genomics Consortium meta-analysis of major depression comprises 688,808 cases and 4,364,225 controls across 29 countries, identifying 697 independent associations at 635 loci¹¹ — and *SLC1A2* is not among them, a point conceded in the limitations section of the positive study itself.⁹

We therefore asked a narrow, tractable question: what does the *SLC1A2*/*SLC1A3* depression literature actually contain, and does the positive finding survive forensic examination using only information already in the public record?

## 2. Methods

### 2.1 Search and eligibility

We searched for case-control association studies reporting genotype counts for *SLC1A2* or *SLC1A3* variants in MDD, across eight search strategies: rs4354668 × MDD; other *SLC1A2* variants × MDD; *SLC1A3* × MDD; post-stroke depression; antidepressant response and treatment resistance; suicidal behaviour; rs4354668 × schizophrenia; and non-European candidate-gene panels. Reference lists of all eligible studies and of two recent reviews were hand-searched. Studies testing within-patient outcomes (treatment response, grey-matter volume, spectroscopic glutamate, cognition) or lacking a control group were excluded from the primary analysis and retained for cross-disorder context. Reporting follows PRISMA 2020¹² and the STREGA extension for genetic association studies;¹³ the flow of studies is given in Figure 3 and Supplementary Table S1.

### 2.2 Data harmonisation

rs4354668 appears in the literature under three encodings (A/C, T/G, −181A>C) reflecting opposite-strand reporting. The alleles are non-palindromic, so harmonisation is unambiguous: A ≡ T, C ≡ G, with C/G the low-activity allele. All counts were harmonised to A/C.

### 2.3 Statistical analysis

Odds ratios were computed under allelic, recessive and dominant models with Woolf standard errors. Pooling used inverse-variance fixed effects and DerSimonian–Laird random effects¹⁴; because DerSimonian–Laird is severely anticonservative at small *k*, we report the Hartung–Knapp *t*-interval¹⁵ as the primary random-effects result, per current guidance for *k* < 5.¹⁶ Hardy–Weinberg equilibrium was assessed by χ² goodness-of-fit; heterozygote deficit was quantified as Wright's *F* with 4,000-draw multinomial bootstrap intervals. Population differentiation used the two-population Wright's *F*ST.¹⁷ Power was computed for the allelic χ² test; false-positive report probability followed Wacholder et al.¹⁸ evaluated at plausible rather than observed effect sizes. Excess significance followed Ioannidis and Trikalinos.¹⁹ Single-variant discrimination was computed analytically from genotype distributions under a log-additive model.

Because cell counts are small, genotypic contrasts were additionally tested by Fisher's exact test with a mid-*p* correction. Bayesian reanalysis used a normal prior on the log odds ratio with standard deviation calibrated to observed common-variant effects in depression, reporting the posterior mean and the Bayes factor against the null.²⁰ Quantitative bias analysis modelled differential genotyping error as per-allele mis-scoring and solved for the rate required to generate the observed case distribution from a true distribution equal to the controls.²¹ Vibration of effects enumerated every odds ratio obtainable from the same dataset across five inheritance models and three strata.²² A Monte Carlo simulation (40,000 replicates) modelled the published record under a range of true effects and suppression probabilities for non-significant findings. Leave-one-out and cumulative analyses are reported for completeness. Analyses used Python 3.13 (numpy 2.4, scipy 1.17); all code and extracted counts are provided (§ Code availability).

### 2.4 Use of a large language model

A large language model (Claude, Anthropic) was used to execute the search strategies listed in §2.1, to extract genotype counts from the published tables of the included studies, to write and run the analysis code described in §2.3, and to prepare the first draft of this manuscript. Every extracted count was checked against the corresponding source publication, and every reported statistic is reproduced by the deposited code. The model is not an author and holds no responsibility for the content. The authors accept full responsibility for the accuracy and integrity of the work.

## 3. Results

### 3.1 The evidence base is two studies

Eight search strategies identified **two** eligible MDD case-control studies of rs4354668 (Table 1, Figure 3), totalling **261 cases and 562 controls**. For *SLC1A3* and MDD diagnosis, no eligible study exists; the nearest neighbours examine suicidal behaviour or schizophrenia. The originating authors state independently that no further rs4354668 MDD studies have been reported, which our search confirms.

For scale, the entire global rs4354668 MDD literature constitutes approximately 1/2,600th of a single current PGC analysis.

**Table 1. Eligible and contextual case-control studies of rs4354668.**

| Study | Phenotype | Ancestry | Cases | Controls | Counts | Result |
|---|---|---|---|---|---|---|
| Rodek 2022⁹ | MDD | European (Poland) | 161 | 462* | Yes | Positive |
| Thaweethee 2019²³ | MDD | SE Asian (Thailand) | 100 | 100 | Yes | Null |
| Merk 2019²⁴ | Schizophrenia | European (Poland) | 393 | 462* | Yes | Null |
| Zhang 2015²⁵ | Schizophrenia | E Asian (China) | ~1,117 | ~2,289 | No | Positive |
| Nagai 2009²⁶ | Schizophrenia | E Asian (Japan) | ~200 | ~200 | No | Null |
| Deng 2004²⁷ | Schizophrenia | E Asian (Japan) | ~100 | ~100 | No | Null |

\* Identical control series (see §3.3).

### 3.2 Pooling produces an uninformative interval

Fixed-effect pooling of the two MDD studies gives an allelic OR of 1.45 (95% CI 1.16–1.81, *p* = 0.0009) — an apparently publishable result. The individual estimates were OR 1.57 (95% CI 1.21–2.02) in the Polish sample and 1.16 (0.75–1.78) in the Thai sample. The Hartung–Knapp interval, the only defensible random-effects summary at *k* = 2,¹⁶ spans **0.23 to 8.67** (Figure 1a; Table 2). The data are compatible with a strong protective effect, no effect, or a strong risk effect.

The two studies are additionally poorly poolable: control C-allele frequencies differ by 28 percentage points (0.413 in Poland, 0.695 in Thailand; two-population *F*ST = 0.080). Under such differentiation, linkage disequilibrium with any true causal variant will differ between samples even if the causal effect is identical. Publication bias could not be assessed; Egger, Begg, trim-and-fill and funnel inspection all require *k* ≥ 10.

**Table 2. Pooled estimates, MDD only.**

| Model | Fixed (IV) | Random (DL) | **Random (Hartung–Knapp)** | Q | I² |
|---|---|---|---|---|---|
| Allelic | 1.45 (1.16–1.81) | 1.42 (1.07–1.88) | **1.42 (0.23–8.67)** | 1.40 | 29% |
| Recessive | 1.88 (1.35–2.63) | 1.80 (1.05–3.07) | **1.80 (0.06–57.8)** | 2.40 | 58% |
| Dominant | 1.37 (0.94–2.00) | 1.37 (0.94–2.00) | **1.37 (0.15–12.7)** | 0.83 | 0% |

### 3.3 The two European studies share a control series

Merk et al.²⁴ (schizophrenia) and Rodek et al.⁹ (MDD) report identical control genotype counts (AA 154, AC 234, CC 74; n = 462), the same recruitment site, the same bioethics approval number, the same published primer pair and the same restriction enzyme. These are one control series, not two. Any synthesis treating them as independent would double-count 462 participants — a duplication invisible from the abstracts and detectable only by comparing raw tables.

The overlap also creates an unusually clean natural experiment.

### 3.4 The finding fails internal replication

Against the same controls, with the same assay, in the same laboratory and city, rs4354668 is null in **2.4-fold more patients** (Table 3, Figure 1a).

**Table 3. Same controls, same assay, two patient series.**

| Group | n | CC (%) | C allele (%) | Allelic OR | Recessive OR | *p* |
|---|---|---|---|---|---|---|
| Controls | 462 | 16.0 | 41.3 | ref | ref | — |
| **Schizophrenia** | **393** | **19.3** | **43.3** | **1.08 (0.89–1.31)** | **1.26 (0.88–1.79)** | **0.20** |
| MDD | 161 | 30.4 | 52.5 | 1.57 (1.21–2.02) | 2.29 (1.51–3.48) | 0.0001 |

This is materially more informative than an external replication failure. The most frequent explanation for spurious candidate-gene results is a poorly matched control group — here, super-selected blood donors 17 years younger on average than cases. That explanation is largely excluded: the same controls behaved as expected against a better-powered patient series. The anomaly localises to the MDD case series.

A direct case-case comparison, which removes the controls from the analysis entirely, confirms this: CC genotype, MDD versus schizophrenia, **OR 1.82 (1.20–2.77), *p* = 0.005**. Two psychiatric inpatient groups from the same region genotyped on the same bench differ from each other more than the schizophrenia group differs from healthy donors. Given the substantial genetic correlation between these disorders and the cross-diagnostic character of common psychiatric risk variants,²⁸ disorder specificity of this magnitude would be a substantial outlier.

### 3.5 Genotyping quality control

All series pass Hardy–Weinberg testing (*p* > 0.05), a weak screen. Wright's *F* is more sensitive to the dominant PCR-RFLP failure mode, incomplete digestion, which miscalls heterozygotes as homozygotes.

The MDD case series shows the largest heterozygote deficit of any series in the dataset — *F* = +0.116 (95% CI −0.039 to +0.273) — exceeding the schizophrenia cases run on the same bench (+0.026) and the shared controls (−0.044); the two Thai series were unremarkable (cases +0.022, controls −0.109) (Figure 1b). The interval crosses zero, so this is **hypothesis-generating rather than probative**. We state it as such deliberately: over-interpreting an interval that includes the null would replicate the error this paper documents. What can be said is that the one series producing a positive result also carries the pattern differential genotyping error would produce, in the direction that would generate a spurious recessive-model association.

### 3.6 The reported sex-specific effect does not survive testing

Rodek et al. reported significance in women (*p* < 0.001) and non-significance in men (*p* = 0.052), and constructed a mechanistic account around it: oestrogen upregulates GLT-1, so a low-expression genotype might attenuate oestrogen-mediated neuroprotection.

That inference compares *p*-values across strata. Comparing the effect estimates gives: women OR 2.59 (1.51–4.44), men OR 2.08 (1.01–4.26), **ratio of ORs 1.25 (0.51–3.06), interaction z = 0.48, *p* = 0.63**. There is no evidence of sex modification; the estimates are nearly identical and differ only in precision (116 women versus 45 men). This is the difference-between-significant-and-non-significant fallacy.²⁹ The mechanism invoked is real biology, deployed to explain an artefact of subgroup size — a route by which candidate-gene literatures acquire spurious mechanistic depth.

### 3.7 Power, multiplicity and credibility

At the effect sizes observed for common variants in depression (OR 1.02–1.05), the positive study had **3.5–5.7% power** — indistinguishable from the type I error rate. Detecting OR 1.05 at 80% power requires 13,574 cases at α = 0.05, or 68,488 at genome-wide significance. False-positive report probability at a generous prior of 0.05 is 0.39 assuming a true OR of 1.20 and 0.91 assuming 1.05 (Table 4).

The reported *p* survives Bonferroni correction across the full 45–105-test analysis grid (adjusted *p* = 0.003–0.008). It does not survive the selection step: the authors used Bayesian information criterion to choose the best-fitting inheritance model, then reported inference from the selected model on the same data without penalty for the selection.

**Table 4. False-positive report probability, positive recessive finding (*p* ≈ 0.001).**

| Prior | True OR 1.05 | True OR 1.20 |
|---|---|---|
| 0.25 | 0.63 | 0.09 |
| 0.05 | 0.91 | 0.39 |
| 0.01 | 0.98 | 0.77 |

### 3.8 The literature contains excess significance

Across the six located studies, two significant findings were reported. Expected numbers given each study's power (Figure 1c):

| Assumed true OR | Expected | Observed | Binomial *p* |
|---|---|---|---|
| 1.02 | 0.23 | 2 | **0.019** |
| 1.05 | 0.41 | 2 | **0.059** |
| 1.10 | 0.93 | 2 | 0.237 |
| 1.20 | 2.19 | 2 | 0.707 |
| 1.45 | 4.29 | 2 | 0.991 |

At the effect sizes modern psychiatric genetics observes, these studies had collective power to expect well under one significant result. Two were published. The excess resolves only if the true OR is assumed to be ≥1.20 — i.e. only by assuming the hypothesis under test, and at a magnitude the PGC data have effectively excluded for common depression variants.

By contrast, the direction of effect is consistent across all four studies reporting one (sign test *p* = 0.125). We note that with *k* = 4 this test cannot reach significance even at 4/4, and that a literature reporting only confirmatory analyses generates directional consistency by construction.

### 3.9 Clinical utility is capped independently of validity

The proposed diagnostic application can be evaluated arithmetically without resolving the association. A single biallelic variant at frequency 0.41 with the published per-allele OR of 1.57 yields **AUC 0.584** (Figure 1d). Chance is 0.500; screening instruments are generally considered unusable below ~0.75; the polygenic score derived from the full genome-wide analysis, aggregating millions of variants, achieves only modest discrimination.³⁰˒¹¹ Taking the recessive OR of 2.38 at face value, genotyping moves lifetime risk from 16% to 39% while missing 68% of future cases.

This is a ceiling imposed by effect size, not precision. No accumulation of further studies alters it.

---

### 3.10 Bayesian reanalysis

A *p*-value does not give the probability that an association is real. Placing a prior on the log odds ratio calibrated to the effect sizes common variants actually show in depression (standard deviation 0.05, so that 95% of true effects lie between OR 0.91 and 1.10), the published estimate of 1.57 shrinks to a posterior of **1.06** (Figure 2a), an 87% reduction. The Bayes factor against the null is **2.0** — "barely worth a bare mention" on the conventional scale — despite an observed *p* below 0.001. Only under a prior permitting effects up to OR 2.7, which no contemporary genome-wide analysis supports for a common variant, does the Bayes factor exceed 50 (Supplementary Table S8).

### 3.11 Genotyping error does not explain the finding

The heterozygote deficit reported in §3.5 suggests differential genotyping error, and we tested that explanation rather than asserting it.

Modelling the error as per-allele mis-scoring, and taking the true case distribution to equal the control distribution, the rate required to generate the observed excess of CC homozygotes is **24.5%**. At that rate the 5% re-genotyping the authors report would have yielded approximately 99% discordance; they report 100% concordance. The one-parameter model also fails to reproduce the observed distribution even at that rate, under-predicting AA by a third (goodness of fit χ² = 4.9, 2 df). At plausible rates of 1–5%, mis-scoring generates apparent recessive odds ratios of 1.04 to 1.20, against an observed 2.29 (Figure 2b; Supplementary Table S9).

**Differential genotyping error is therefore insufficient**, and we do not advance it as the explanation. The anomaly localises to the case series by §3.4, but its source is more plausibly case selection, cryptic population substructure, or chance.

### 3.12 Vibration of effects

Fifteen defensible analyses of the same dataset — five inheritance models across three strata — yield odds ratios from **0.64 to 2.59**, a four-fold range, with 9 of 15 reaching *p* < 0.05 unadjusted (Figure 2c; Supplementary Table S10). The reported result is not the largest available. This spread is what a reader requires in order to interpret any single selected estimate, and it is why the inheritance model must be fixed in advance rather than chosen by fit.

### 3.13 The published record is consistent with a null plus selective reporting

Modelling publication bias as suppression of non-significant findings, we simulated the fraction of published studies reaching significance across a range of true effects (Figure 2d; Supplementary Table S11). Under a true null with complete reporting, two or more significant results among six published studies arises in 0.9% of simulated literatures. Under a true null with 80–90% of negative findings unpublished, it arises in 14% to 35%.

The observed record is therefore consistent with no effect plus selective reporting, and also with a small true effect of OR 1.05–1.10 under complete reporting. It is **not** consistent with the published effect size: at OR 1.20 or above, substantially more significant findings would be expected than were observed. A literature of six small studies can exclude a large effect; it cannot distinguish a small one from none.

## 4. Discussion

Three independent lines of evidence, derived from different data, converge on the same conclusion. The published data cannot support the claim (§3.1–3.2, 3.7). The claim fails internal replication within the originating laboratory (§3.4). And the literature's own pattern of published results is statistically anomalous (§3.8). Notably, the argument strengthens as it proceeds: the first line is an absence of evidence, the second is contradictory evidence, and the third is positive evidence of reporting bias.

**The functional-validation trap.** rs4354668's appeal rests on genuine functional validation.⁸ We suggest this has been systematically over-weighted. Demonstrating that a variant alters promoter activity establishes that it does *something*; it does not privilege any particular downstream phenotype. rs4354668 has now been tested against progressing stroke,⁸ reward dependence,³¹ schizophrenia (four times),²⁴˒²⁵˒²⁶˒²⁷ bipolar illness course,³² hippocampal grey matter,³³ working memory,³⁴ Parkinson's disease³⁵ and depression.⁹˒²³ Functional plausibility multiplied across many outcomes generates positives by chance, and functionality is precisely what licenses the multiplication.

We were unable to locate any report of rs4354668 as a brain expression quantitative trait locus for *SLC1A2* in GTEx, PsychENCODE, or published multi-ancestry brain eQTL meta-analyses. This is an absence of found evidence rather than a demonstrated null, and we flag it as a priority check rather than a finding. But it identifies a gap at the base of the literature: the mechanistic premise of every study cited here derives from a reporter construct in cell culture, not from *in vivo* expression in brain. If the promoter effect does not survive native chromatin context, the mechanistic account offered by these papers is unsupported in the tissue that matters.

**A conceptual tension worth resolving.** This literature routinely pairs "glutamatergic hypofunction" with "impaired clearance" in the same paragraph. These pull in opposite directions: reduced EAAT2 increases, not decreases, extracellular glutamate. The reconciliation lies in compartmentalisation. Impaired clearance drives spillover onto extrasynaptic NMDA receptors, whose activation produces CREB shut-off and BDNF suppression, while degrading the signal-to-noise of phasic synaptic transmission.³⁶ The result is simultaneously excess tonic glutamate and hypofunctional phasic signalling — and it accounts for the otherwise paradoxical efficacy of an NMDA antagonist. This is the defensible form of the astrocytic hypothesis, and it is interrogable by magnetic resonance spectroscopy and positron emission tomography rather than by genotyping.

**Two anticipated objections.** First, that genome-wide studies rely on minimal phenotyping, whose specificity for MDD is low, so a genome-wide null need not apply to rigorously diagnosed samples.³⁷ This argument has real force in general, and Merk et al. make a version of it explicitly about their homogeneous Silesian sample.²⁴ It is nonetheless answerable, because the current PGC release separates clinical- and interview-ascertained cohorts into their own summary statistics; that subset matches the DSM-5 inpatient ascertainment of the study examined here and should be interrogated as a co-primary analysis. Second, that the focus on rs4354668 may be misplaced. In the only study to apply multiple-testing correction across a panel of *SLC1A2* and *SLC1A3* variants, the sole SNP to survive was rs12294045; the nominal *SLC1A3* signal at rs2731880 did not.³⁸ If rs12294045 carries the real signal, two decades of work have been directed at the wrong position in the gene, and it has not been followed up. Those authors also concede in their limitations that neither gene appears among the 287 loci of the current schizophrenia GWAS³⁹ — an admission parallel to the one made by the authors of the positive depression study.⁹

**On the source study.** Our critique is methodological and implies nothing about the integrity of the investigators. Rodek and colleagues published complete genotype counts — which is what made every analysis in this paper possible — reported the contradictory Thai result rather than omitting it, reused their control series transparently across two publications, and stated in their own limitations that *SLC1A2* lies outside GWAS-indexed depression loci. That is more candour than much of this literature offers. The five-model grid, sex stratification, haplotype analysis and clinical-variable ANOVA are all standard practice. That is the point: standard practice produced 5-HTTLPR,⁴⁰ and it is still in use.

**Generalisability of the audit.** The six checks applied here require only published tables:

1. Compare raw control counts across papers from the same group to detect sample overlap.
2. Compute Wright's *F* per series; HWE testing alone is insufficient.
3. Search for same-laboratory studies of related disorders sharing the control series.
4. Test case-case where two patient series share genotyping.
5. Test subgroup claims as interactions, never as differences in significance.
6. Apply the excess-significance test across the located literature.

None requires raw data, ethics approval or new recruitment. Applied prospectively at peer review, this sequence would have flagged the present claim before publication.

**Implications for the field.** Border et al.¹⁰ argued in 2019 that candidate-gene hypotheses for depression should be abandoned. The literature examined here was published in 2022 and continues to be cited as supporting a diagnostic biomarker. Formal recommendations to abandon an approach evidently do not, by themselves, stop its outputs entering the record. Journal-level requirements — mandatory *a priori* power calculations against contemporary effect-size expectations, mandatory lookup of the focal variant in the relevant consortium summary statistics, and mandatory pre-registration of the genetic model — would be more effective, and cost a referee minutes.

### 4.1 Limitations

First, we could not extract genotype counts from three East Asian schizophrenia studies; these contribute to the vote-count and excess-significance analyses but not to pooled estimates, and their sample sizes are approximate. Second, the heterozygote-deficit observation is suggestive only. Third, we did not query PGC summary statistics directly; our statements about genome-wide evidence rest on the published loci lists and on the source authors' own concessions. Fourth, excess-significance testing is itself low-powered at *k* = 6 and sensitive to the assumed true effect, which we address by reporting across a range. Fifth, our search may have missed non-indexed or non-English studies, though the resulting bias would most likely be toward unpublished nulls, strengthening rather than weakening the conclusion.

### 4.2 Conclusion

*SLC1A2* rs4354668 is not a supported susceptibility locus for major depressive disorder, and could not serve as a diagnostic marker even if it were. The underlying astrocytic glutamate hypothesis remains well supported by post-mortem, pharmacological and cell-type-level genomic evidence, and is best advanced through expression, imaging and Mendelian randomisation designs. The definitive genetic test — variant and gene-level interrogation of existing consortium data comprising 2,600-fold more cases than the entire candidate literature — costs an afternoon and has not been done.

---

## Declarations

**Funding.** No funding was received for conducting this study.

**Competing interests.** A.K. and T.K. hold invention disclosures and pending patents in two technology areas: macromolecular structures for neurological indications, and responsive biomaterials for microenvironment-targeted delivery. None of these disclosures or applications covers *SLC1A2*, *SLC1A3*, the variants examined here, or the methods applied in this work. A.K. is a Novo Nordisk AI Fellow at MIT (from November 2026). T.K. holds an Anthropic AI for Science award supporting work other than this study; the large language model used in this work (§2.4) is a product of Anthropic, disclosed here for completeness. Neither the fellowship nor the award provided financial, computational, or operational support for this study.

**Ethics approval.** Not applicable. This work is a secondary analysis of aggregate genotype counts published in the cited studies. No individual-level data, human participants or animals were involved.

**Consent to participate.** Not applicable.

**Consent for publication.** Not applicable.

**Data availability.** All genotype counts analysed here are reproduced verbatim from the published tables of the cited studies and are tabulated in full in Supplementary Table S2. No unpublished or individual-level data were used.

**Code availability.** All analysis code, with source genotype counts inline, is provided as Supplementary Software and will be deposited at Zenodo on acceptance [DOI to be assigned].

**Author contributions.** A.K. and T.K. contributed equally to this work. Conceptualisation: A.K., T.K. Methodology: A.K., T.K. Software: A.K., T.K. Formal analysis: A.K., T.K. Data curation: A.K., T.K. Writing — original draft: A.K., T.K. Writing — review and editing: A.K., T.K. Both authors read and approved the final manuscript.

**Acknowledgements.** We thank the authors of the primary studies for publishing complete genotype tables, without which this analysis would not have been possible.

---

## References

1. Rajkowska G, Stockmeier CA. Astrocyte pathology in major depressive disorder: insights from human postmortem brain tissue. *Curr Drug Targets* 2013;14:1225–36.
2. Choudary PV, Molnar M, Evans SJ, et al. Altered cortical glutamatergic and GABAergic signal transmission with glial involvement in depression. *Proc Natl Acad Sci USA* 2005;102:15653–8.
3. Medina A, Burke S, Thompson RC, et al. Glutamate transporters: a key piece in the glutamate puzzle of major depressive disorder. *J Psychiatr Res* 2013;47:1150–6.
4. Chandley MJ, Szebeni K, Szebeni A, et al. Gene expression deficits in pontine locus coeruleus astrocytes in men with major depressive disorder. *J Psychiatry Neurosci* 2013;38:276–84.
5. Bechtholt-Gompf AJ, Walther HV, Adams MA, et al. Blockade of astrocytic glutamate uptake in rats induces signs of anhedonia and impaired spatial memory. *Neuropsychopharmacology* 2010;35:2049–59.
6. John CS, Smith KL, Van't Veer A, et al. Blockade of astrocytic glutamate uptake in the prefrontal cortex induces anhedonia. *Neuropsychopharmacology* 2012;37:2467–75.
7. Als TD, Kurki MI, Grove J, et al. Depression pathophysiology, risk prediction of recurrence and comorbid psychiatric disorders using genome-wide analyses. *Nat Med* 2023;29:1832–44.
8. Mallolas J, Hurtado O, Castellanos M, et al. A polymorphism in the EAAT2 promoter is associated with higher glutamate concentrations and higher frequency of progressing stroke. *J Exp Med* 2006;203:711–7.
9. Rodek P, Kowalczyk M, Kowalski J, et al. Association study of the SLC1A2 (rs4354668), SLC6A9 (rs2486001), and SLC6A5 (rs2000959) polymorphisms in major depressive disorder. *J Clin Med* 2022;11:5914.
10. Border R, Johnson EC, Evans LM, et al. No support for historical candidate gene or candidate gene-by-interaction hypotheses for major depression across multiple large samples. *Am J Psychiatry* 2019;176:376–87.
11. Adams MJ, et al.; Major Depressive Disorder Working Group of the Psychiatric Genomics Consortium. Trans-ancestry genome-wide study of depression identifies 697 associations implicating cell types and pharmacotherapies. *Cell* 2025;188:640–652.e9.
12. Page MJ, McKenzie JE, Bossuyt PM, et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. *BMJ* 2021;372:n71.
13. Little J, Higgins JPT, Ioannidis JPA, et al. STrengthening the REporting of Genetic Association studies (STREGA): an extension of the STROBE statement. *PLoS Med* 2009;6:e1000022.
14. DerSimonian R, Laird N. Meta-analysis in clinical trials. *Control Clin Trials* 1986;7:177–88.
15. Hartung J, Knapp G. A refined method for the meta-analysis of controlled clinical trials with binary outcome. *Stat Med* 2001;20:3875–89.
16. IntHout J, Ioannidis JPA, Borm GF. The Hartung–Knapp–Sidik–Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian–Laird method. *BMC Med Res Methodol* 2014;14:25.
17. Weir BS, Cockerham CC. Estimating F-statistics for the analysis of population structure. *Evolution* 1984;38:1358–70.
18. Wacholder S, Chanock S, Garcia-Closas M, et al. Assessing the probability that a positive report is false: an approach for molecular epidemiology studies. *J Natl Cancer Inst* 2004;96:434–42.
19. Ioannidis JPA, Trikalinos TA. An exploratory test for an excess of significant findings. *Clin Trials* 2007;4:245–53.
20. Goodman SN. Toward evidence-based medical statistics. 2: The Bayes factor. *Ann Intern Med* 1999;130:1005–13.
21. Lash TL, Fox MP, MacLehose RF, et al. Good practices for quantitative bias analysis. *Int J Epidemiol* 2014;43:1969–85.
22. Patel CJ, Burford B, Ioannidis JPA. Assessment of vibration of effects due to model specification can demonstrate the instability of observational associations. *J Clin Epidemiol* 2015;68:1046–58.
23. Thaweethee B, Suttajit S, Thanoi S, et al. Association of SLC1A2 and SLC17A7 polymorphisms with major depressive disorder in a Thai population. *Asian Biomed* 2019;12:131–8.
24. Merk W, Kucia K, Mędrala T, et al. Association study of the excitatory amino acid transporter 2 (EAAT2) and glycine transporter 1 (GlyT1) gene polymorphism with schizophrenia in a Polish population. *Neuropsychiatr Dis Treat* 2019;15:989–1000.
25. Zhang B, Guan F, Chen G, et al. Common variants in SLC1A2 and schizophrenia: association and cognitive function in patients with schizophrenia and healthy individuals. *Schizophr Res* 2015;169:128–34.
26. Nagai Y, Ohnuma T, Karibe J, et al. No genetic association between the SLC1A2 gene and Japanese patients with schizophrenia. *Neurosci Lett* 2009;463:223–7.
27. Deng X, Shibata H, Ninomiya H, et al. Association study of polymorphisms in the excitatory amino acid transporter 2 gene (SLC1A2) with schizophrenia. *BMC Psychiatry* 2004;4:21.
28. Cross-Disorder Group of the Psychiatric Genomics Consortium. Genetic relationship between five psychiatric disorders estimated from genome-wide SNPs. *Nat Genet* 2013;45:984–94.
29. Gelman A, Stern H. The difference between "significant" and "not significant" is not itself statistically significant. *Am Stat* 2006;60:328–31.
30. Wray NR, Ripke S, Mattheisen M, et al. Genome-wide association analyses identify 44 risk variants and refine the genetic architecture of major depression. *Nat Genet* 2018;50:668–81.
31. Matsumoto Y, Suzuki A, Ishii G, et al. The −181 A/C polymorphism in the excitatory amino acid transporter-2 gene promoter affects the personality trait of reward dependence in healthy subjects. *Neurosci Lett* 2007;427:99–102.
32. Dallaspezia S, Poletti S, Lorenzi C, et al. Influence of an interaction between lithium salts and a functional polymorphism in SLC1A2 on the history of illness in bipolar disorder. *Mol Diagn Ther* 2012;16:303–9.
33. Poletti S, Locatelli C, Radaelli D, et al. Effect of early stress on hippocampal gray matter is influenced by a functional polymorphism in EAAT2 in bipolar disorder. *Prog Neuropsychopharmacol Biol Psychiatry* 2014;51:146–52.
34. Spangaro M, Bosia M, Zanoletti A, et al. Cognitive dysfunction and glutamate reuptake: effect of EAAT2 polymorphism in schizophrenia. *Neurosci Lett* 2012;522:151–5.
35. Appenzeller S, Schulte C, Thier S, et al. No association between polymorphisms in the glutamate transporter SLC1A2 and Parkinson's disease. *Mov Disord* 2013;28:1305–6.
36. Hardingham GE, Bading H. Synaptic versus extrasynaptic NMDA receptor signalling: implications for neurodegenerative disorders. *Nat Rev Neurosci* 2010;11:682–96.
37. Flint J. The genetic basis of major depressive disorder. *Mol Psychiatry* 2023;28:2254–65.
38. Wang L, Ma T, Qiao D, et al. Polymorphism of rs12294045 in EAAT2 gene is potentially associated with schizophrenia in Chinese Han population. *BMC Psychiatry* 2022;22:171.
39. Trubetskoy V, Pardiñas AF, Qi T, et al. Mapping genomic loci implicates genes and synaptic biology in schizophrenia. *Nature* 2022;604:502–8.
40. Culverhouse RC, Saccone NL, Horton AC, et al. Collaborative meta-analysis finds no evidence of a strong interaction between stress and 5-HTTLPR genotype contributing to the development of depression. *Mol Psychiatry* 2018;23:133–42.


## Figure legends

**Figure 1. Forensic summary.** **(a)** Allelic odds ratios for rs4354668. The fixed-effect pooled estimate excludes unity; the Hartung–Knapp interval appropriate at *k* = 2 spans an order of magnitude in each direction. Asterisks denote the shared control series. **(b)** Wright's *F* by series with bootstrap intervals; the single positive series carries the largest heterozygote deficit. **(c)** Observed versus expected significant findings across six studies. Shading marks the range of effect sizes consistent with contemporary genome-wide evidence. **(d)** Analytic ceiling on single-variant discrimination.

**Figure 2. Extended analyses.** **(a)** Posterior odds ratio under normal priors of increasing width, with the Bayes factor against the null above each estimate. Under the prior calibrated to common-variant effects in depression (SD 0.05) the published estimate shrinks by 87%. **(b)** Apparent recessive odds ratio generated from a true null by per-allele mis-scoring. Shading marks plausible genotyping error rates; reproducing the observed estimate requires 24.5%. **(c)** Every odds ratio obtainable from the same dataset across five inheritance models and three strata; open circles are not significant at *p* < 0.05. **(d)** Probability that at least two of six published studies reach significance, by true effect size and by *q*, the probability that a non-significant study is published.

**Figure 3. PRISMA 2020 flow of studies.** Counts are reported from full-text assessment onward, and every assessed record is named in Supplementary Table S1. Four reports excluded from the primary analysis are retained for the cross-disorder comparison and the excess-significance analysis.


---
---

# ⚠ PRE-SUBMISSION NOTES — NOT PART OF THE MANUSCRIPT

## The one thing that must be added before submission

**Run the PGC lookup and put it in the Results.** As written, the paper says "the definitive test costs an afternoon and has not been done" — and then does not do it. The first referee will ask why. Adding a Results subsection reporting rs4354668's β, SE and *p* in PGC MDD 2025, PGC3 SCZ and PGC BD, plus a two-one-sided-test equivalence result establishing |OR| < 1.02, converts the paper from *"this claim is unsupported"* to *"this claim is refuted."* That single addition is likely worth a journal tier. The same applies to the eQTL lookup (§4).

Both require downloading public summary statistics — no permissions, no new data.

## Realistic target journals

| Tier | Venue | Fit |
|---|---|---|
| Strong | *American Journal of Psychiatry* | Published Border 2019; natural home for a psychiatric-genetics correction |
| Strong | *Molecular Psychiatry* — Expert Review or Brief Comm. | Receptive to meta-research on candidate genes |
| Good | *Biological Psychiatry* / *Biol Psychiatry GOS* | Fits the methodological brief |
| Good | *Translational Psychiatry* | Fast, open access |
| Alternative | *Royal Society Open Science* / *Meta-Psychology* | If framed primarily as the reusable audit framework |

*Nature*/*Science*/*Cell* are not realistic: this is a negative-result and meta-research paper, not new biology. The framing that maximises reach is the **generalisable audit framework** (§4 "Generalisability"), with rs4354668 as the worked case — consider promoting that into the title and abstract if targeting a methods-oriented venue.

## Other pre-submission items

- **Right of reply.** Consider notifying the Rodek and Merk groups before submission and offering co-authorship on a joint re-analysis. Several journals now expect this for papers critiquing identified work, and a jointly authored correction lands better than an adversarial one.
- **Soften §3.5 further if a referee objects.** The heterozygote-deficit interval crosses zero. It is currently framed as hypothesis-generating; if pressed, it can be moved to Supplementary without weakening the paper, since §3.4 carries the argument.
- **PRISMA 2020 flow diagram** and **PROSPERO registration** should be added — reviewers will require both for anything framed as a systematic review.
- **Verify the PGC 2025 case count** (688,808 vs 685,808 appears in preprint and final versions) against the version of record before submission.
- **Confirm approximate sample sizes** for Nagai 2009 and Deng 2004 from the primary papers; these are currently estimates and are flagged as such in Limitations.
