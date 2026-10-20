# Provenance

## Parent release

- Parent repository: `peng777wang/cross-cancer-direction-v1`
- Parent release: `v1.3.0`, published 2026-10-16
- Concept DOI: `10.5281/zenodo.23141927`
- Exact v1.3.0 DOI: `10.5281/zenodo.23416675`
- GitHub release commit: the tag `v1.3.0` on `main`

The parent release already contains the 273-locus resource, the replication
panel resource, simulation results, colocalization results and the interactive
locus map. This repository does not duplicate those large resources.

## Post-v1.3.0 artifacts

The following analyses were performed on 2026-10-19 after the v1.3.0 release:

1. `f07_susie_p0.R` attempted fine-mapping at 42 reversed loci and produced 84
   locus-arm fits. The reversed arms returned no credible set; the concordant
   arms returned singleton credible sets and `susieR` reported inconsistent
   summary statistics and LD.
2. The f07 result was rejected as evidence and documented rather than discarded.
3. `f08_region_zscan.R` replaced the fit with a ±500 kb regional z scan that
   does not use an LD matrix.
4. The 42-row scan was formatted as Supplementary Table 12 for the manuscript.

## Verification performed for this package

- f07 output rows: 84
- f07 lead-concordance rows: 42
- f08 output rows: 42
- Supplementary Table 12 rows: 42
- Supplementary Table 12 records above `|z| = 5.45`: 0
- Supplementary Table 12 records marked `N`: 42

File-level hashes are in `metadata/SHA256SUMS.txt`.

## Independence from the GPB submission package

This repository is an analysis addendum. It contains no manuscript PDF, cover
letter, reviewer response, Editorial Manager export, or other GPB submission
file. Those files remain in the local submission workspace and are not part of
this public code and data package.
