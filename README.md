# Effects of impaired glycemic homeostasis on specific regional changes in brain structure: evidence from a Mendelian randomization study

R code accompanying a bidirectional, two-sample Mendelian randomization (MR) study of four glycemic traits and 78 MRI-derived brain structural phenotypes.

## Data and setup

- **Glycemic traits:** fasting glucose (FG), fasting insulin (FI), glycated hemoglobin (HbA1c), and 2-hour post-challenge glucose (2hGlu), from [MAGIC](https://magicinvestigators.org/downloads/index.html).
- **Brain phenotypes:** 35 cortical surface-area traits, 35 cortical-thickness traits, seven subcortical volumes, and intracranial volume, from [ENIGMA](https://enigma.ini.usc.edu/research/download-enigma-gwas-results/). Some ENIGMA summary statistics require an approved access request.

Obtain the GWAS releases used in the study and prepare the local input files expected by `Code.R`, including `FG.tsv`, `FI.tsv`, `HbA1c.tsv`, `2hGlu.tsv`, `SNP_Traitsrow.csv`, and the brain-phenotype files named in the script. The summary statistics are not included here. Install the R packages `MendelR` and `FastTraitR` and their dependencies, then run the script from the directory containing the prepared inputs:

```bash
Rscript Code.R
```

## Analysis workflow

1. **Screen and prepare instruments.** The study uses European-ancestry GWAS data, independent exposure-associated SNPs, allele harmonization, instrument-strength and Steiger checks, and review of potential pleiotropic associations. The relevant thresholds and screening decisions are detailed in the manuscript and its supplementary tables.
2. **Run forward MR.** `mr_common_batch()` evaluates each glycemic trait against the 78 brain phenotypes in three groups: cortical surface area, cortical thickness, and volumes. The script requests local LD clumping and excludes SNPs strongly associated with the outcome (`p2 = 5e-8`).
3. **Check pleiotropy and re-estimate selected results.** `mr_presso_batch()` analyzes the three forward-result files. The subsequent `mr_common()` calls repeat selected exposure–outcome analyses after the SNP exclusions listed in the script.
4. **Run reverse MR.** The final three `mr_common_batch()` calls use brain phenotypes as exposures and glycemic traits as outcomes, with a relaxed instrument threshold (`p1 = 5e-6`).
5. **Interpret multiple tests and sensitivity analyses.** The manuscript evaluates the inverse-variance weighted estimates across 312 forward comparisons using Benjamini–Hochberg FDR correction, alongside alternative MR estimators, statistical power, and additional sensitivity analyses.

`Code.R` contains the batch MR, MR-PRESSO, selected SNP-exclusion reruns, and reverse-MR calls. Manual GWAS Catalog screening, FDR adjustment, power calculations, and PCMR are reported in the manuscript but are not separate executable steps in this file. Consult the manuscript and supplementary tables for the full methods and results.

## Citation

Please cite the associated article if you use this code or its results.

