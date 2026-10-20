# Cross-cancer direction audit addendum

This repository is a standalone addendum to the v1.3.0
`cross-cancer-direction-v1` analysis release. It records the abandoned SuSiE
fine-mapping attempt, the LD-independent regional z scan that replaced it, and
the source table used for Supplementary Table 12.

It is intentionally separate from the Genomics, Proteomics & Bioinformatics
(GPB) manuscript submission package. No journal manuscript, cover letter, or
Editorial Manager material is included here.

## Result

- The SuSiE analysis was not used in the manuscript.
- Across the 42 reversed loci, no reversed-arm regional variant with a sign
  opposing the MTAG lead reached the genome-wide threshold of `|z| = 5.45`.
- The largest opposing `|z|` was 4.53. All 42 loci are recorded as `N` in
  Supplementary Table 12.
- The regional scan requires no LD matrix and is therefore not exposed to the
  reference-panel mismatch that invalidated the SuSiE run.

## Repository contents

| Path | Purpose |
|---|---|
| `docs/ANALYSIS_LOG_f07_susie_2026-10-19.md` | Six-part technical record: attempt, outcome, diagnosis, decision, replacement and artifacts |
| `docs/PROVENANCE.md` | Relationship to the v1.3.0 release and file-level provenance |
| `scripts/f07_susie_p0.R` | SuSiE fit retained for auditability; its result is not used in the manuscript |
| `scripts/f08_region_zscan.R` | LD-independent replacement analysis |
| `source_data/f07_p0_credsets.csv` | 84 fitted locus-arm records from f07 |
| `source_data/f07_lead_concordance.csv` | Lead-variant concordance audit |
| `source_data/f08_region_zscan.csv` | 42-row regional scan result |
| `source_data/Supplementary_Table_12_region_zscan.csv` | Source table used for Supplementary Table 12 |
| `logs/f07_run.log` | Run log, including the `susieR` LD-consistency warning |
| `metadata/SHA256SUMS.txt` | SHA-256 checksums for the distributed files |
| `metadata/manifest.json` | Machine-readable package inventory and key checks |

## Reproduce the regional scan

Set `CGS_ROOT` to a directory containing the inputs used by the parent release:

```text
loci273.txt
pleiotropy_resource.csv
f07_p0_credsets.csv
reg_<locus>_z.tsv
```

Then run:

```bash
Rscript scripts/f08_region_zscan.R
```

The script writes `f08_region_zscan.csv`. The packaged result is the output used
to build Supplementary Table 12.

For the archived SuSiE audit, set `CGS_ROOT` and optionally `CGS_R_LIB`, then:

```bash
Rscript scripts/f07_susie_p0.R
```

This rerun is provided for provenance only. Its degenerate fits should not be
interpreted as credible-set evidence.

## Parent analysis release

- Concept DOI: `10.5281/zenodo.23141927` (resolves to v1.3.0)
- Exact v1.3.0 version DOI: `10.5281/zenodo.23416675`
- Parent repository: `https://github.com/peng777wang/cross-cancer-direction-v1`

## Data and attribution

The package contains derived results from publicly available cancer GWAS
summary statistics. No individual-level data are included. Input datasets and
accessions are documented in Supplementary Table 1 of the associated
manuscript.

If you use this addendum, cite it together with the parent analysis release
above.

## Licence

MIT. See `LICENSE`.
