#!/global/home/hpc5173/miniconda3/envs/main/bin/Rscript

# -----------------------------------------------------------------------------
# Title: Boruta-based feature importance for CRISPR dependency by cancer type
# Author: Michael Vermeulen
# Purpose: For each input gene, run Boruta random-forest feature selection to
#          rank mutation features (damaging + hotspot flags) associated with
#          CRISPR gene dependency (DepMap 23Q2), stratified by cancer type.
#
# Inputs:
#   - args[1]: path to a text file with one gene symbol per line.
#              Only genes found in the CRISPRGeneDependency matrix are kept.
#
# Data dependencies (update paths as needed):
#   - DepMap/CCLE 23Q2 CSVs:
#       * OmicsSomaticMutationsMatrixDamaging.csv
#       * OmicsSomaticMutationsMatrixHotspot.csv
#       * CRISPRGeneDependency.csv (used as response)
#       * Model.csv (for OncotreePrimaryDisease annotations)
#     (exp/CN/fusions/PRISM loaded below for completeness but not used here.)
#
# Outputs (per cancer type and gene):
#   - /global/home/hpc5173/projects/SL_network/output/RF/MUT_DEPENDENCY_23Q2/<TYPE>/<GENE>_getImpExtraGini.RDS
#   - ..._getImpExtraRaw.RDS
#   - ..._getImpRfZ.RDS
#     Each file stores Boruta attStats() with an added `normImp` column and a
#     `cancer_type` label equal to the output directory basename.
#
# Notes:
#   - Cancer-type filter: computed from available models; types with >8 lines are kept.
#     (The original inline comment said ">=19" but the code uses >8; kept as >8 here.)
#   - Parallelization uses pbmclapply (forking). Works on Linux/macOS; on Windows
#     pbmcapply may fall back to sequential execution. Tune mc.cores below.
#   - Idempotency: existence of the *_getImpExtraGini.RDS acts as a sentinel to
#     skip recomputation for a given (TYPE, GENE). Delete to re-run.
#   - Stochasticity: Boruta/RandomForest are stochastic. Consider `set.seed()`
#     before Boruta calls for reproducibility.
# -----------------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
args[1] -> ind  # path to gene list (one symbol per line)

# ---- Libraries ----
library(tidyverse)
library(Boruta)
library(data.table)
library(stringr)
library(pbmcapply)

# ---- Read input gene list ----
ind %>% as.character() -> ind
# one column, no header
data.table::fread(file = ind, header = FALSE) -> n
# head(n, 20)  # optional peek

# ---- Load DepMap / PRISM 23Q2 resources ----
# Core matrices
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsSomaticMutationsMatrixDamaging.csv", nThread = 8) -> dam
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsSomaticMutationsMatrixHotspot.csv", nThread = 8) -> hot
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsExpressionProteinCodingGenesTPMLogp1.csv", nThread = 8) -> exp       # loaded, not used
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/CRISPRGeneEffect.csv", nThread = 8) -> effect                              # loaded, not used
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsCNGene.csv", nThread = 8) -> CN                                      # loaded, not used
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Model.csv", nThread = 8) -> model
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsFusionFiltered.csv", nThread = 8) -> fusion                           # loaded, not used
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/CRISPRGeneDependency.csv", nThread = 8) -> dep
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv", nThread = 8) -> PRISM
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv", nThread = 8) -> PRISM_info  # loaded, not used

# ---- Normalize column names (strip trailing " (…)") ----
stringr::str_split(string = colnames(dam),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(dam)
stringr::str_split(string = colnames(hot),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(hot)
stringr::str_split(string = colnames(exp),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(exp)
stringr::str_split(string = colnames(effect), pattern = " \\(", simplify = TRUE)[,1] -> colnames(effect)
stringr::str_split(string = colnames(CN),     pattern = " \\(", simplify = TRUE)[,1] -> colnames(CN)
stringr::str_split(string = colnames(dep),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(dep)
stringr::str_split(string = colnames(PRISM),  pattern = " \\(", simplify = TRUE)[,1] -> drug_names  # parsed for reference only

# ---- Response: dependency (not effect) ----
effect <- dep  # prefer CRISPRGeneDependency as response

# ---- Restrict to genes present in dependency matrix ----
n$V1[n$V1 %in% colnames(effect)] -> g

# ---- Prepare PRISM matrix (parsed but unused) ----
tibble::column_to_rownames(PRISM, 'V1') -> PRISM
PRISM %>% t() %>% as.data.frame() -> PRISM

# ---- Build mutation feature matrix: damaging + hotspot ----
tibble::column_to_rownames(dam, 'V1') -> dam
tibble::column_to_rownames(hot, 'V1') -> hot

colnames(dam) <- paste0(colnames(dam), "_DAM")
colnames(hot) <- paste0(colnames(hot), "_HOT")

dplyr::left_join(
  x = dam %>% tibble::rownames_to_column("ID"),
  y = hot %>% tibble::rownames_to_column("ID"),
  by = "ID"
) %>% tibble::column_to_rownames("ID") -> tmp

# Replace NAs with 0 (no mutation)
tmp[is.na(tmp)] <- 0

# For joins later
tibble::rownames_to_column(PRISM, 'V1') -> PRISM  # not used beyond this
names(effect)[1] <- "V1"  # make model ID explicit for joins

# ---- Select cancer types with sufficient sample size (>8 lines) ----
# Use TP53 dependency to anchor model IDs (any gene column would work)
effect %>% dplyr::select(V1, "TP53") -> tm
names(tm)[2] <- "eff"

dplyr::left_join(y = tmp %>% tibble::rownames_to_column("V1"), x = tm, by = 'V1') -> df

dplyr::left_join(
  x = df,
  y = model %>% dplyr::select(1, "OncotreePrimaryDisease"),
  by = c("V1" = "ModelID")
) -> df

# Keep complete cases, then count per cancer type
df[complete.cases(df)] -> df

df %>%
  dplyr::group_by(OncotreePrimaryDisease) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  dplyr::arrange(dplyr::desc(count)) -> samples

samples[samples$count > 8, ]$OncotreePrimaryDisease -> types  # vector of cancer types to analyze

# ---- Main loop over cancer types ----
lapply(X = types, FUN = function(cancer_types){

  message(cancer_types)

  # Build output directory name like "NonSmallCellLungCancer" (no spaces/hyphens)
  gsub(x = cancer_types, pattern = " ", replacement = "") -> pathname
  gsub(x = pathname,     pattern = "-", replacement = "") -> pathname
  pathname <- paste0("/global/home/hpc5173/projects/SL_network/output/RF/MUT_DEPENDENCY_23Q2/", pathname)
  if (!dir.exists(pathname)) { dir.create(pathname, recursive = TRUE) }

  # Parallel over genes
  pbmclapply(X = g, mc.cores = 8, FUN = function(gene){

    message(gene)

    # Idempotency sentinel
    if (file.exists(paste0(pathname, "/", gene, "_getImpExtraGini.RDS"))) { return(NULL) }

    # ---- Build design frame for this (type, gene) ----
    # Response y: dependency for the current gene
    effect %>% dplyr::select(V1, dplyr::all_of(gene)) -> dd
    names(dd)[2] <- "eff"

    # Predictors X: mutation features (damaging + hotspot)
    dplyr::left_join(y = tmp %>% tibble::rownames_to_column("V1"), x = dd, by = 'V1') -> df4

    # Add cancer type
    dplyr::left_join(
      y = df4,
      x = model %>% dplyr::select("ModelID", "OncotreePrimaryDisease"),
      by = c("ModelID" = "V1")
    ) -> df4

    # Keep complete cases and filter to the current type
    df4[complete.cases(df4)] -> df4
    df4 %>% dplyr::filter(OncotreePrimaryDisease %in% cancer_types) -> df4

    # Drop non-informative numeric columns (all-zero predictors)
    df4 <- df4 %>% dplyr::select(c(1:3), where(~ is.numeric(.) && sum(.) != 0))

    # ---- Boruta with three importance measures ----
    # Helper: train, convert to attStats, add normImp + label, save RDS
    run_and_save <- function(getImpFun, suffix){
      fit <- Boruta(
        eff ~ .,
        data    = df4 %>% dplyr::filter() %>% dplyr::select(-1, -2),  # drop IDs & type
        pValue  = 0.01,
        getImp  = getImpFun,
        maxRuns = 500
      )
      out <- attStats(fit)
      out %>% dplyr::mutate(normImp = medianImp / sum(medianImp)) -> out
      out$cancer_type <- basename(pathname)
      saveRDS(out, file = paste0(pathname, "/", gene, suffix, ".RDS"))
    }

    run_and_save(getImpExtraGini, "_getImpExtraGini")
    run_and_save(getImpExtraRaw,  "_getImpExtraRaw")
    run_and_save(getImpRfZ,       "_getImpRfZ")

    return(NULL)
  }) -> n  # list of NULLs; not used

  return(NULL)

}) -> nn  # list of NULLs; not used

