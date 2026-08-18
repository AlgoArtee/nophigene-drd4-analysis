# Changelog

## Unreleased

- Refactored genetic variant interpretation to decode sample-level VCF `GT`, zygosity, allele dosage, and genotype-call QC before phenotype inference.
- Added configurable confidence heuristics for FILTER, QUAL, DP, AD balance, GQ, PL/GP support, strand bias, and orientation-bias evidence.
- Updated predictive theses and reports to separate site-level `REF -> ALT` definitions from sample genotype state, with conservative HERC2/OCA2 eye-colour uncertainty wording.
