# Genotype-Aware Variant Analysis

## Why GT is authoritative

VCF `REF` and `ALT` describe the alleles that can exist at a genomic site. They do not, by themselves, describe the sample's genotype. The sample-level `GT` field is the authoritative field for genotype state.

Examples:

- `REF=A`, `ALT=G`, `GT=0/0` decodes to `A/A`.
- `REF=A`, `ALT=G`, `GT=0/1` decodes to `A/G`.
- `REF=A`, `ALT=G`, `GT=1/1` decodes to `G/G`.
- `REF=A`, `ALT=G,T`, `GT=1/2` decodes to `G/T`.

The analysis therefore treats `REF -> ALT` as a site definition only. Phenotype-oriented text is based on decoded `GT`, zygosity, ALT dosage, and call quality.

## Why INFO AF/AC/AN are not substitutes

VCF `INFO` fields such as `AC`, `AF`, and `AN` are site-level summary fields. In cohort VCFs they summarize many samples. In single-sample VCFs they may summarize that file's call set, but they still do not replace the sample genotype.

Sample-specific evidence belongs in the FORMAT/sample columns, including `GT`, `AD`, `DP`, `GQ`, `PL`, `GP`, `AF`, `SB`, `F1R2`, and `F2R1`. This app uses `GT` to decode genotype and uses those supporting fields to explain confidence.

## Dosage, QC, And Uncertainty

Trait interpretation now distinguishes homozygous reference, heterozygous, homozygous alternate, and multiallelic genotypes. A heterozygous genotype is not collapsed into the same evidence state as a homozygous alternate genotype.

The call-confidence layer flags missing GT, non-PASS FILTER, low depth, low GQ, weak PL/GP support, heterozygous allelic imbalance, DP versus AD inconsistency, and strand/orientation-bias evidence when available. These QC flags lower confidence in the genotype call; they do not rewrite the genotype.

For polygenic traits such as HERC2/OCA2 eye colour, predictions stay probabilistic. For example, rs12913832 is a major contributor, but `A/G` carries less certainty than `G/G`, and brown or hazel outcomes can remain plausible when only a small SNP subset is available.
