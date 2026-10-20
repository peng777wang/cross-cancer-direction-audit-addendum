# Documented technical issue: SuSiE fine-mapping of the reversed loci

## Attempted

Fine-map each of the 42 reversed loci in its MTAG-significant arms with
`susie_rss` (two-stage signal number: L = 3, escalated to L = 10 on
saturation), using the 1000G-EUR reference panel already used at 15q25.
Script: `f07_susie_p0.R`. Output: `f07_p0_credsets.csv`.

## Outcome

84 fits (42 loci x 2 arms). The reversed arms (gastric and small-cell lung)
returned **no credible set at all** in 42 of 42 fits (top PIP 0.005-0.06). The
concordant arms (prostate) returned 3-10 credible sets, but **every credible set
was a single variant with PIP = 1.000**, and 41 of 84 fits escalated to L = 10
only to return ten singleton sets. `susieR` emitted its standard diagnostic:

```
Please check consistency between summary statistics and LD matrix.
```

## Diagnosis

Two separate problems, only one of which is expected.

* **Expected**: the reversed arm has no own-GWAS signal, so there is nothing to
  localize. This reproduces the 0-of-42 result already reported and adds no
  information.
* **Not usable**: the prostate-arm output is the classic signature of a
  summary-statistic / LD mismatch - every variant appears independent and each
  receives PIP = 1. This is the **same class of silent mismatch** the manuscript
  already guards against at 15q25, where the Methods state that "subsetting the
  table before indexing the matrix is a silent mismatch". Here the exposure is
  the reference panel rather than an index, but the failure mode is identical:
  `susie_rss` does not raise an error, it returns a degenerate but plausible
  looking result.

## Decision

The fine-mapping run is **not used in the manuscript**. Comparing the lead
variant across arms - the diagnostic originally planned - is uninformative by
construction, because one arm has no credible set at all; the resulting 0-of-42
concordance would be an artifact of an empty fit rather than evidence about
localization.

## Replacement

A regional z scan (`f08_region_zscan.R`) tests the same hypothesis without LD:
for each reversed arm, across a +/-500 kb window around the lead variant, it records the largest
variant with a sign opposing the MTAG lead. Results are reported in
Supplementary Table 12 and in Results: at all 42 loci the largest opposing
|z| was 4.53, below the genome-wide threshold of 5.45 and within the range
expected by chance for regions of this size.

## Artifacts

* `f07_susie_p0.R`, `f07_run.log`, `f07_p0_credsets.csv`, `f07_lead_concordance.csv`
* `f08_region_zscan.R`, `f08_region_zscan.csv`
* `Supplementary_Table_12_region_zscan.csv`
