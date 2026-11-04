# ============================================================================
# Study title:
# Effects of impaired glycemic homeostasis on specific regional changes in
# brain structure: evidence from a Mendelian randomization study
#
# GWAS summary statistics: original data sources
# Exposures: MAGIC (Meta-Analyses of Glucose and Insulin-related traits
# Consortium): fasting glucose (FG), 2-hour post-challenge glucose (2hGlu),
# glycated hemoglobin (HbA1c), and fasting insulin (FI).
# Download portal: https://magicinvestigators.org/downloads/index.html
#
# Outcomes: ENIGMA Consortium: cortical surface area, cortical thickness,
# subcortical brain volumes, and intracranial volume.
# Data access and download portal:
# https://enigma.ini.usc.edu/research/download-enigma-gwas-results/
#
# IMPORTANT: Select the GWAS release matching the datasets used in the study.
# Some ENIGMA summary statistics require a data-access request. Input names
# below refer to locally prepared files rather than original download names.
# ============================================================================

# Run this script from the directory containing the prepared GWAS input files.
# setwd("path/to/prepared/GWAS/files")
getwd()
library(MendelR)
library(FastTraitR)
res =look_trait(file_name="SNP_Traitsrow.csv", pval=1e-5)

# Forward MR: glycemic traits and cortical surface area
mr_common_batch(c("FG.tsv","2hGlu.tsv","HbA1c.tsv","FI.tsv"), c("SurfArea","SA-bankssts","SA-caudalanteriorcingulate","SA-caudalmiddlefrontal",
                                                                                  "SA-cuneus","SA-entorhinal","SA-frontalpole","SA-fusiform","SA-inferiorparietal",
                                                                                  "SA-inferiortemporal","SA-insula","SA-isthmuscingulate","SA-lateraloccipital",
                                                                                  "SA-lateralorbitofrontal","SA-lingual","SA-medialorbitofrontal","SA-middletemporal",
                                                                                  "SA-paracentral","SA-parahippocampal","SA-parsopercularis","SA-parsorbitalis",
                                                                                  "SA-parstriangularis","SA-pericalcarine","SA-postcentral","SA-posteriorcingulate",
                                                                                  "SA-precentral","SA-precuneus","SA-rostralanteriorcingulate","SA-rostralmiddlefrontal",
                                                                                  "SA-superiorfrontal","SA-superiorparietal","SA-superiortemporal","SA-supramarginal",
                                                                                  "SA-temporalpole","SA-transversetemporal"),local_clump=T,p2 = 5e-08)

mr_presso_batch("01.glycemic_35SA.csv") # Input file for batch MR-PRESSO analysis
# Forward MR: glycemic traits and cortical thickness
mr_common_batch(c("FG.tsv","2hGlu.tsv","HbA1c.tsv","FI.tsv"), c("Thickness","TH-bankssts","TH-caudalanteriorcingulate","TH-caudalmiddlefrontal","TH-cuneus",
                                                                                  "TH-entorhinal","TH-frontalpole","TH-fusiform","TH-inferiorparietal",
                                                                                  "TH-inferiortemporal","TH-insula","TH-isthmuscingulate","TH-lateraloccipital",
                                                                                  "TH-lateralorbitofrontal","TH-lingual","TH-medialorbitofrontal","TH-middletemporal",
                                                                                  "TH-paracentral","TH-parahippocampal","TH-parsopercularis","TH-parsorbitalis",
                                                                                  "TH-parstriangularis","TH-pericalcarine","TH-postcentral","TH-posteriorcingulate",
                                                                                  "TH-precentral","TH-precuneus","TH-rostralanteriorcingulate","TH-rostralmiddlefrontal",
                                                                                  "TH-superiorfrontal","TH-superiorparietal","TH-superiortemporal",
                                                                                  "TH-supramarginal","TH-temporalpole","TH-transversetemporal"),local_clump=T,p2 = 5e-08)
mr_presso_batch("01.glycemic_35TH.csv") # Input file for batch MR-PRESSO analysis
# Forward MR: glycemic traits and subcortical volumes
mr_common_batch(c("FG.tsv","2hGlu.tsv","HbA1c.tsv","FI.tsv"), c("Nucleus accumbens_1042","Amygdala_1043","Caudate_1044","Pallidum_1046","Putamen_1047",
                                                                                  "Thalamus_1048","Hippocampus_1045","Intracranial_1041"),local_clump=T,p2 = 5e-08)
mr_presso_batch("01.glycemic_volume.csv") # Input file for batch MR-PRESSO analysis

################################### Outlier exclusion: cortical surface area
a <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs12692738") 
mr_common("2hGlu.tsv", "SA-entorhinal",local_clump=T,rm_snps = a)  ### One outcome-specific outlier excluded: "rs12692738"

a <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs12692738") 
mr_common("2hGlu.tsv", "SA-inferiortemporal",local_clump=T,rm_snps = a)  ### One outcome-specific outlier excluded: "rs12692738"

a <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs117643180","rs2126259","rs878521") 
mr_common("2hGlu.tsv", "SA-middletemporal",local_clump=T,rm_snps = a)  ### Three outcome-specific outliers excluded: "rs117643180","rs2126259","rs878521"

a <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs17271305") 
mr_common("2hGlu.tsv", "SA-rostralanteriorcingulate",local_clump=T,rm_snps = a)  ### One outcome-specific outlier excluded: "rs17271305"

a <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs10838693") 
mr_common("FG.tsv", "SA-paracentral",local_clump=T,rm_snps = a)  ### One outcome-specific outlier excluded: "rs10838693"

a <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs7012637") 
mr_common("FG.tsv", "SA-superiortemporal",local_clump=T,rm_snps = a)  ### One outcome-specific outlier excluded: "rs7012637"

a <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs10774624") 
mr_common("HbA1c.tsv", "SA-bankssts",local_clump=T,rm_snps = a)  ### One outcome-specific outlier excluded: "rs10774624"
#"FG.tsv","2hGlu.tsv","HbA1c.tsv","FI.tsv"
################################### Outlier exclusion: cortical thickness
b <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs12692738") 
mr_common("2hGlu.tsv", "TH-posteriorcingulate",local_clump=T,rm_snps = b)  ### One outcome-specific outlier excluded: "rs12692738"

b <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs550057") 
mr_common("2hGlu.tsv", "TH-superiorfrontal",local_clump=T,rm_snps = b)  ### One outcome-specific outlier excluded: "rs550057"

b <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs550057") 
mr_common("2hGlu.tsv", "TH-superiortemporal",local_clump=T,rm_snps = b)  ### One outcome-specific outlier excluded: "rs550057"

b <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs7012637") 
mr_common("FG.tsv", "TH-fusiform",local_clump=T,rm_snps = b)  ### One outcome-specific outlier excluded: "rs7012637"

b <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs10838524","rs10974438") 
mr_common("FG.tsv", "TH-insula",local_clump=T,rm_snps = b)  ### Two outcome-specific outliers excluded: "rs10838524","rs10974438"

b <- c("rs1260326","rs13389219","rs2001846","rs4148646","rs1800562","rs10865959",
       "rs7012814","rs13089972","rs1800437","rs11610045","rs2649999","rs7903146",
       "rs1820176") 
mr_common("FG.tsv", "TH-isthmuscingulate",local_clump=T,rm_snps = b)  ### One outcome-specific outlier excluded: "rs1820176"

### Reverse-direction Mendelian randomization analysis
# Reverse MR: cortical surface area and glycemic traits
mr_common_batch(c("SurfArea","SA-bankssts","SA-caudalanteriorcingulate","SA-caudalmiddlefrontal",
                                                                "SA-cuneus","SA-entorhinal","SA-frontalpole","SA-fusiform","SA-inferiorparietal",
                                                                "SA-inferiortemporal","SA-insula","SA-isthmuscingulate","SA-lateraloccipital",
                                                                "SA-lateralorbitofrontal","SA-lingual","SA-medialorbitofrontal","SA-middletemporal",
                                                                "SA-paracentral","SA-parahippocampal","SA-parsopercularis","SA-parsorbitalis",
                                                                "SA-parstriangularis","SA-pericalcarine","SA-postcentral","SA-posteriorcingulate",
                                                                "SA-precentral","SA-precuneus","SA-rostralanteriorcingulate","SA-rostralmiddlefrontal",
                                                                "SA-superiorfrontal","SA-superiorparietal","SA-superiortemporal","SA-supramarginal",
                                                                "SA-temporalpole","SA-transversetemporal"),c("FG.tsv","2hGlu.tsv","HbA1c.tsv","FI.tsv"), p1=5e-06)
# Reverse MR: cortical thickness and glycemic traits
mr_common_batch(c("Thickness","TH-bankssts","TH-caudalanteriorcingulate","TH-caudalmiddlefrontal","TH-cuneus",
                                                                "TH-entorhinal","TH-frontalpole","TH-fusiform","TH-inferiorparietal",
                                                                "TH-inferiortemporal","TH-insula","TH-isthmuscingulate","TH-lateraloccipital",
                                                                "TH-lateralorbitofrontal","TH-lingual","TH-medialorbitofrontal","TH-middletemporal",
                                                                "TH-paracentral","TH-parahippocampal","TH-parsopercularis","TH-parsorbitalis",
                                                                "TH-parstriangularis","TH-pericalcarine","TH-postcentral","TH-posteriorcingulate",
                                                                "TH-precentral","TH-precuneus","TH-rostralanteriorcingulate","TH-rostralmiddlefrontal",
                                                                "TH-superiorfrontal","TH-superiorparietal","TH-superiortemporal",
                                                                "TH-supramarginal","TH-temporalpole","TH-transversetemporal"),c("FG.tsv","2hGlu.tsv","HbA1c.tsv","FI.tsv"),  p1=5e-06)
# Reverse MR: subcortical volumes and glycemic traits
mr_common_batch( c("Nucleus accumbens_1042","Amygdala_1043","Caudate_1044","Pallidum_1046","Putamen_1047",
                                                                "Thalamus_1048","Hippocampus_1045","Intracranial_1041"),c("FG.tsv","2hGlu.tsv","HbA1c.tsv","FI.tsv"), p1=5e-06)
