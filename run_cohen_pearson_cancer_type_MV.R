#!/global/home/hpc5173/miniconda3/envs/main/bin/Rscript

# -----------------------------------------------------------------------------
# Title: Cohen's d and Pearson correlation for dependency vs mutation by tissue
# Author: Michael Vermeulen
# Purpose: For (source,target,tissue) triples nominated by previous Boruta runs,
#          compute (i) Pearson correlation between mutation status (0/1) of the
#          source gene and CRISPR dependency of the target gene, and (ii)
#          Cohen's d effect size comparing target dependency between mutated vs
#          wild-type groups, stratified by tissue (OncotreePrimaryDisease) or PAN.
#
# Inputs:
#   - args[1]: path to a text file with one source gene symbol per line. These
#              are used to subset rows from the wide `test` table (from Boruta).
#
# Data dependencies (update paths as needed):
#   - DepMap/CCLE 23Q2 CSVs: damaging & hotspot mutation matrices, CRISPR
#     gene dependency/effect matrices, Model.csv (for tissue labels), etc.
#   - Precomputed Boruta summary RDS (wide format per source/target/tissue):
#       /global/home/hpc5173/projects/SL_network/output/RF/MUT_DEPENDENCY_23Q2/
#         Tissue_important_DEP_all.RDS
#     Expected to contain columns: source, target, cancer_type, test, medianImp.
#
# Outputs (per source gene):
#   - /global/home/hpc5173/projects/SL_network/output/COHENS/DEP/
#       Tissue_important_DEP_<SOURCE>_EFFECT.RDS
#     An RDS data.frame joining Boruta stats for (source,target,tissue) with
#     computed pearson r/p and Cohen's d (+ 95% CI), plus group counts.
#
# Notes:
#   - This script treats the source mutation feature as binary: >0 -> 1 else 0.
#   - It computes Pearson correlation between the binary source indicator and
#     the continuous target dependency; and Cohen's d between groups.
#   - The code expects the `effsize` package for `cohen.d()`. Install/load it.
#       install.packages("effsize"); library(effsize)
#   - IMPORTANT: `dep` is currently loaded from CRISPRGeneEffect.csv. If you
#     intended to use CRISPRGeneDependency.csv (as elsewhere), update that path.
#   - Special-case: the tissue label "BLymphoblasticLeukemia" is mapped to
#     "BLymphoblasticLeukemia/Lymphoma" to match Model.csv annotations.
#   - Parallelization over (source,target,tissue) uses pbmclapply (forking).
# -----------------------------------------------------------------------------

## add/data
library(tidyverse)
library(Boruta)      
library(data.table)
library(stringr)
library(pbmcapply)
library(magrittr)
library(effsize)   

# ---- CLI args ----
args = commandArgs(trailingOnly=TRUE)
args[1] -> ind   # path to a text file with one source gene per line

# Read the input gene list (no header)
ind %>% as.character() -> ind
data.table::fread(file = ind, header = FALSE) -> n

# ---- Load DepMap / PRISM resources ----
# Mutation feature matrices
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsSomaticMutationsMatrixDamaging.csv", nThread = 8) -> dam
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsSomaticMutationsMatrixHotspot.csv",  nThread = 8) -> hot
# Optional additional layers (loaded but not used downstream here)
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsExpressionProteinCodingGenesTPMLogp1.csv", nThread = 8) -> exp
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/CRISPRGeneEffect.csv",                                 nThread = 8) -> effect
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsCNGene.csv",                                       nThread = 8) -> CN
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Model.csv",                                             nThread = 8) -> model
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsFusionFiltered.csv",                               nThread = 8) -> fusion
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/CRISPRGeneEffect.csv",                                 nThread = 8) -> dep
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv", nThread = 8) -> PRISM
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv",  nThread = 8) -> PRISM_info

# ---- Normalize column names (strip trailing " (…)") to align model IDs ----
stringr::str_split(string = colnames(dam),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(dam)
stringr::str_split(string = colnames(hot),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(hot)
stringr::str_split(string = colnames(exp),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(exp)
stringr::str_split(string = colnames(effect), pattern = " \\(", simplify = TRUE)[,1] -> colnames(effect)
stringr::str_split(string = colnames(CN),     pattern = " \\(", simplify = TRUE)[,1] -> colnames(CN)
stringr::str_split(string = colnames(dep),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(dep)
stringr::str_split(string = colnames(PRISM),  pattern = " \\(", simplify = TRUE)[,1] -> drug_names  # kept for reference

# ---- Build mutation feature matrix and join dependency/effect ----
# Index by model ID
tibble::column_to_rownames(dam, 'V1') -> dam
tibble::column_to_rownames(hot, 'V1') -> hot

# Suffix columns to distinguish damaging vs hotspot flags
colnames(dam) <- paste0(colnames(dam), "_DAM")
colnames(hot) <- paste0(colnames(hot), "_HOT")

# Join damaging + hotspot by model ID -> `tmp` (rowname is model ID)
dplyr::left_join(
  x = dam %>% tibble::rownames_to_column("ID"),
  y = hot %>% tibble::rownames_to_column("ID"),
  by = "ID"
) %>% tibble::column_to_rownames("ID") -> tmp

# Join dependency/effect matrix by model ID (adds target dependency columns)
dplyr::left_join(
  x = tmp %>% tibble::rownames_to_column("ID"),
  y = dep,
  by = c("ID" = "ModelID")
) -> tmp

# ---- Load Boruta summary and reshape to wide ----
# Expect columns: source, target, cancer_type, test, medianImp
readRDS(file = paste0("/global/home/hpc5173/projects/SL_network/output/RF/MUT_DEPENDENCY_23Q2/Tissue_important_DEP_all.RDS")) -> df

# Drop duplicate (cancer_type, source, target, test) rows, then cast test->columns
dplyr::distinct(df, cancer_type, source, target, test, .keep_all = TRUE) -> df
reshape2::dcast(
  data   = df,
  formula = source + target + cancer_type ~ test,
  value.var = "medianImp"
) -> test

# ---- Annotate with tissue labels and clean names ----
data <- dplyr::left_join(
  x = model %>% dplyr::select(ModelID, OncotreePrimaryDisease),
  y = tmp,
  by = c("ModelID" = "ID")
)
# Normalize tissue label formatting to match later comparisons
gsub(x = data$OncotreePrimaryDisease, pattern = " ", replacement = "") -> data$OncotreePrimaryDisease
gsub(x = data$OncotreePrimaryDisease, pattern = "-", replacement = "") -> data$OncotreePrimaryDisease

# ---- Main statistic function -------------------------------------------------
call <- function(source, target, tissue, data){
  # Handle special-case label used in Model.csv
  if (tissue == "BLymphoblasticLeukemia") {
    tissue <- "BLymphoblasticLeukemia/Lymphoma"
  }

  # Subset matrix to the tissue or all (PAN); keep only IDs + source + target
  if (tissue == "PAN") {
    data %>% dplyr::filter() %>% dplyr::select(ModelID, dplyr::all_of(source), dplyr::all_of(target)) %>% as.data.frame() -> w
  } else {
    data %>% dplyr::filter(OncotreePrimaryDisease == tissue) %>% dplyr::select(ModelID, dplyr::all_of(source), dplyr::all_of(target)) %>% as.data.frame() -> w
  }

  # Drop NAs in either column; binarize source mutation feature (>0 -> 1 else 0)
  w %>% tidyr::drop_na(dplyr::all_of(target), dplyr::all_of(source)) -> w
  ifelse(w[,2] > 0, 1, 0) -> w[,2]

  # Pearson correlation: binary source vs continuous target dependency
  cor.test(
    x = w %>% dplyr::pull(dplyr::all_of(source)),
    y = w %>% dplyr::pull(dplyr::all_of(target)),
    method = "pearson"
  ) -> o

  # Cohen's d: target dependency difference (mutated vs wild-type)
  effsize::cohen.d(
    w[w[,2] == 1, ] %>% dplyr::pull(dplyr::all_of(target)),
    w[w[,2] == 0, ] %>% dplyr::pull(dplyr::all_of(target))
  ) -> oo

  # Group sizes
  w[w[,2] == 1, ] %>% nrow() -> n
  w[w[,2] == 0, ] %>% nrow() -> nn

  # Return a tidy one-row data.frame with stats and counts
  data.frame(
    source       = source,
    target       = target,
    tissue       = tissue,
    pear_c       = o$estimate,
    pear_p       = o$p.value,
    cohens       = oo$estimate,
    cohens_95conf = paste(oo$conf.int, collapse = ":"),
    mut_n        = n,
    n            = n + nn
  ) %>% return()
}
# -----------------------------------------------------------------------------

# Optional: if you had an older RDS path/shape, these show how to rebuild `test`:
# readRDS(file = paste0("/global/home/hpc5173/projects/SL_network/output/RF/Tissue_important_DEP_all.RDS")) -> df
# reshape2::dcast(data = df, formula = source + target + cancer_type  ~ test, value.var = "medianImp") -> test

# Keep rows with all three importance measures present; sort by ExtraGini desc
test %>% tidyr::drop_na(getImpExtraGini, getImpExtraRaw, getImpRfZ) %>% dplyr::arrange(-getImpExtraGini) -> test

# ---- Iterate over input source genes and compute stats per (target,tissue) ----
lapply(X = n$V1 %>% unique(), FUN = function(y){

  # Filter Boruta table to the current source gene
  test %>% dplyr::filter(source == dplyr::all_of(y)) %>% as.data.frame() -> q
  message(y)

  # Parallel over candidate rows for this source gene
  pbmclapply(X = 1:nrow(q), mc.cores = 4, FUN = function(p){
    q[p, ] -> qq

    # Compute stats for this (source,target,tissue)
    call(
      source = qq %>% dplyr::pull(source),
      target = qq %>% dplyr::pull(target),
      tissue = qq %>% dplyr::pull(cancer_type),
      data   = data
    ) -> o

    # Bind Boruta columns with the computed stats (cols 4:9 from `o`)
    cbind(qq, o %>% dplyr::select(4:9)) %>% as.data.frame() %>% return()
  }) %>% data.table::rbindlist() -> df

  # Save per-source RDS
  out_path <- paste0("/global/home/hpc5173/projects/SL_network/output/COHENS/DEP/Tissue_important_DEP_", y, "_EFFECT.RDS")
  message(out_path)
  saveRDS(df, out_path)
  return(NULL)

})

