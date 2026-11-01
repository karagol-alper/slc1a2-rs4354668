# Provenance of the data files

No value in this directory was generated, imputed, estimated or simulated.
Every count is transcribed from a published table.

## genotype_counts.csv

| study_id | Source |
|---|---|
| `rodek22` | Rodek P, Kowalczyk M, Kowalski J, et al. *J Clin Med* 2022;11:5914, **Table 2** |
| `merk19` | Merk W, Kucia K, Mędrala T, et al. *Neuropsychiatr Dis Treat* 2019;15:989–1000, **Table 1** |
| `thaweethee19` | Thaweethee B, Suttajit S, Thanoi S, et al. *Asian Biomed* 2019;12:131–8, **Table 2** |

Counts are stored harmonised to A/C. Thaweethee 2019 publishes T/G; the
mapping is `T -> A`, `G -> C`. The alleles are not palindromic, so this is
unambiguous. The `published_coding` column records the original encoding so
each row can be audited against the source paper.

`eaat.data.load_genotypes` raises if any row's genotypes do not sum to its
declared `n`.

### Shared control series

`rodek22` and `merk19` report the identical control row (AA 154, AC 234,
CC 74; n = 462). Both name the Regional Blood Donation Centre, Katowice, the
same bioethics approval (KNW/0022/KB1/34/14), the same primer pair and the
same *Bcn*I digest. This is one control series, recorded in
`eaat.data.SHARED_CONTROLS`. It is not a transcription error, and the test
suite asserts the two rows stay identical.

## study_ledger.csv

The six located case-control tests of rs4354668 across all psychiatric
phenotypes.

| study_id | Source | Counts published |
|---|---|---|
| `rodek22` | *J Clin Med* 2022;11:5914 | yes |
| `thaweethee19` | *Asian Biomed* 2019;12:131–8 | yes |
| `merk19` | *Neuropsychiatr Dis Treat* 2019;15:989–1000 | yes |
| `zhang15` | *Schizophr Res* 2015;169:128–34 | no |
| `nagai09` | *Neurosci Lett* 2009;463:223–7 | no |
| `deng04` | *BMC Psychiatry* 2004;4:21 | no |

Rows with `n_approximate = 1` publish no genotype counts. Their sample sizes
are estimates taken from the text, and they enter only the vote-count and
excess-significance analyses — never a pooled effect estimate. `control_freq_C`
for those rows is a population-level estimate, not a reported figure.

`significant` records whether the primary rs4354668 comparison was reported as
significant by the original authors, not whether it would survive correction.
