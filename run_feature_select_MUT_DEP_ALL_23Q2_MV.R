#!/global/home/hpc5173/miniconda3/envs/main/bin/Rscript

# -----------------------------------------------------------------------------
# Title: Boruta-based feature importance for CRISPR gene dependency (DepMap 23Q2)
# Author: Michael Vermeulen
# Purpose: For each input gene, run Boruta random-forest feature selection to
#          rank genomic features (damaging & hotspot mutations) associated with
#          CRISPR gene dependency across cancer types with sufficient samples.
#
# Inputs:
#   - args[1]: path to a text file with one gene symbol per line.
#              Only genes present in the dependency matrix are analyzed.
#
# Data dependencies (update paths as needed):
#   - DepMap/CCLE 23Q2 CSVs for mutations (damaging/hotspot), expression (TPM),
#     copy number, CRISPR dependency/effect, model annotations, fusions.
#   - PRISM Repurposing matrices (loaded; column names parsed; not otherwise used).
#
# Outputs (per gene):
#   - <outdir>/<GENE>_getImpExtraGini.RDS
#   - <outdir>/<GENE>_getImpExtraRaw.RDS
#   - <outdir>/<GENE>_getImpRfZ.RDS
#     Each RDS contains Boruta attStats() with an added column `normImp` and the
#     `cancer_type` label (used here only for output directory naming).
#
# Notes:
#   - Parallelization uses pbmclapply (forking). Works on Linux/macOS; on Windows
#     pbmcapply falls back to sequential apply. Adjust mc.cores if needed.
#   - The existence of the *_getImpExtraGini.RDS file is used as a sentinel to
#     skip re-running all three models. Delete it to recompute for a gene.
#   - "PAN" is used as a label for pan-cancer runs (not as a filter). Actual
#     filtering is by cancer types with >8 cell lines (derived below).
#   - Reproducibility: Boruta is stochastic. Consider setting `set.seed()` before
#     the Boruta calls if you want fixed results across runs.
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
# Expecting a single column, no header
data.table::fread(file = ind, header = FALSE) -> n
# head(n, 20)  # quick peek if needed

# ---- Load DepMap / PRISM 23Q2 resources ----
# Large files: prefer fread with threads
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsSomaticMutationsMatrixDamaging.csv", nThread = 8) -> dam
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsSomaticMutationsMatrixHotspot.csv", nThread = 8) -> hot
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsExpressionProteinCodingGenesTPMLogp1.csv", nThread = 8) -> exp
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/CRISPRGeneEffect.csv", nThread = 8) -> effect
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsCNGene.csv", nThread = 8) -> CN
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Model.csv", nThread = 8) -> model
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsFusionFiltered.csv", nThread = 8) -> fusion
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/CRISPRGeneDependency.csv", nThread = 8) -> dep
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv", nThread = 8) -> PRISM
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv", nThread = 8) -> PRISM_info
fread(input = "/global/home/hpc5173/projects/PRISM/data/MSI.csv") -> MSI

# ---- Normalize column names (strip trailing " (…)") for consistency ----
stringr::str_split(string = colnames(dam), pattern = " \\(", simplify = TRUE)[,1] -> colnames(dam)
stringr::str_split(string = colnames(hot), pattern = " \\(", simplify = TRUE)[,1] -> colnames(hot)
stringr::str_split(string = colnames(exp), pattern = " \\(", simplify = TRUE)[,1] -> colnames(exp)
stringr::str_split(string = colnames(effect), pattern = " \\(", simplify = TRUE)[,1] -> colnames(effect)
stringr::str_split(string = colnames(CN), pattern = " \\(", simplify = TRUE)[,1] -> colnames(CN)
stringr::str_split(string = colnames(dep), pattern = " \\(", simplify = TRUE)[,1] -> colnames(dep)
stringr::str_split(string = colnames(PRISM), pattern = " \\(", simplify = TRUE)[,1] -> drug_names  # stored for reference

# ---- Limit to genes present in the dependency matrix ----
# n$V1: gene symbols from the input file
n$V1[n$V1 %in% colnames(effect)] -> g

# ---- Prepare PRISM (parsed but not used downstream) ----
# Kept for compatibility; safe to remove if not needed.
tibble::column_to_rownames(PRISM, 'V1') -> PRISM
PRISM %>% t() %>% as.data.frame() -> PRISM

# ---- Build mutation feature matrix: damaging + hotspot ----
tibble::column_to_rownames(dam, 'V1') -> dam
tibble::column_to_rownames(hot, 'V1') -> hot

colnames(dam) <- paste0(colnames(dam), "_DAM")
colnames(hot) <- paste0(colnames(hot), "_HOT")

# Join on model ID (rowname moved to column `ID`), fill missing with 0 (no mutation)
dplyr::left_join(x = dam %>% tibble::rownames_to_column("ID"),
                 y = hot %>% tibble::rownames_to_column("ID"), by = "ID") %>%
  tibble::column_to_rownames("ID") -> tmp

tmp[is.na(tmp)] <- 0

# For later joins, make sure first column of `effect` is named V1 (model ID)
tibble::rownames_to_column(PRISM, 'V1') -> PRISM
names(effect)[1] <- "V1"

# ---- Identify cancer types with sufficient sample size (>8 lines) ----
# Use TP53 dependency as a quick way to retrieve model IDs (any column would work)
effect %>% dplyr::select(V1, "TP53") -> tm
names(tm)[2] <- "eff"

# Join mutation features to model IDs and annotate with Oncotree primary
# Keep complete rows, then count per cancer type
dplyr::left_join(y = tmp %>% tibble::rownames_to_column("V1"), x = tm, by = 'V1') -> df

dplyr::left_join(x = df,
                 y = model %>% dplyr::select(1, "OncotreePrimaryDisease"),
                 by = c("V1" = "ModelID")) -> df

df[complete.cases(df)] -> df

df %>%
  dplyr::group_by(OncotreePrimaryDisease) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  dplyr::arrange(dplyr::desc(count)) -> samples

samples[samples$count > 8, ]$OncotreePrimaryDisease -> ty  # cancer types to keep

# ---- Run label for pan-cancer (used only for output dir naming) ----
types <- "PAN"

lapply(X = types, FUN = function(cancer_types){

  message(cancer_types)
  gsub(x = cancer_types, pattern = " ", replacement = "") -> pathname
  gsub(x = pathname, pattern = "-", replacement = "") -> pathname
  pathname <- paste0("/global/home/hpc5173/projects/SL_network/output/RF/MUT_DEPENDENCY_23Q2_JAN/", pathname)
  if (!dir.exists(pathname)) { dir.create(pathname, recursive = TRUE) }

  pbmclapply(X = g, mc.cores = 8, FUN = function(gene){

    message(gene)

    # Use the ExtraGini output as a sentinel to skip recomputation
    if (file.exists(paste0(pathname, "/", gene, "_getImpExtraGini.RDS"))) {
      return(NULL)
    }

    # ---- Build design frame for this gene ----
    # Response: CRISPR dependency for `gene` (column renamed to `eff`)
    effect %>% dplyr::select(V1, dplyr::all_of(gene)) -> dd
    names(dd)[2] <- "eff"

    # Predictors: damaging/hotspot mutation flags (joined by model ID)
    dplyr::left_join(y = tmp %>% tibble::rownames_to_column("V1"), x = dd, by = 'V1') -> df4

    # Cancer-type annotation
    dplyr::left_join(y = df4,
                     x = model %>% dplyr::select("ModelID", "OncotreePrimaryDisease"),
                     by = c("ModelID" = "V1")) -> df4

    # Keep complete rows and restrict to cancer types with >8 samples
    df4[complete.cases(df4)] -> df4
    df4 %>% dplyr::filter(OncotreePrimaryDisease %in% ty) -> df4

    # Optional: add MSI, CN, expression, fusions as predictors (commented here)
    # dplyr::left_join(x = df4 , y = MSI, by = c("ModelID" = "DepMap_ID")) -> df4

    # Remove non-informative numeric columns (all zeros)
    df4 <- df4 %>% dplyr::select(c(1:3), where(~ is.numeric(.) && sum(.) != 0))

    # ---- Boruta with three importance measures ----
    # 1) getImpExtraGini
    Boruta(
      eff ~ .,
      data   = df4 %>% dplyr::filter() %>% dplyr::select(-1, -2),  # drop ModelID, Oncotree
      pValue = 0.01,
      getImp = getImpExtraGini,
      maxRuns = 500
    ) -> output

    attStats(output) -> output

    output %>% dplyr::mutate(normImp = medianImp / sum(medianImp)) -> output
    output$cancer_type <- basename(pathname)

    saveRDS(output, file = paste0(pathname, "/", gene, "_getImpExtraGini.RDS"))

    # 2) getImpExtraRaw
    Boruta(
      eff ~ .,
      data   = df4 %>% dplyr::filter() %>% dplyr::select(-1, -2),
      pValue = 0.01,
      getImp = getImpExtraRaw,
      maxRuns = 500
    ) -> output

    attStats(output) -> output

    output %>% dplyr::mutate(normImp = medianImp / sum(medianImp)) -> output
    output$cancer_type <- basename(pathname)

    saveRDS(output, file = paste0(pathname, "/", gene, "_getImpExtraRaw.RDS"))

    # 3) getImpRfZ
    Boruta(
      eff ~ .,
      data   = df4 %>% dplyr::filter() %>% dplyr::select(-1, -2),
      pValue = 0.01,
      getImp = getImpRfZ,
      maxRuns = 500
    ) -> output

    attStats(output) -> output

    output %>% dplyr::mutate(normImp = medianImp / sum(medianImp)) -> output
    output$cancer_type <- basename(pathname)

    saveRDS(output, file = paste0(pathname, "/", gene, "_getImpRfZ.RDS"))

    return(NULL)
  }) -> n  # `n` collects the pbmclapply results (NULLs); not used

  return(NULL)

}) -> nn  # outer lapply result; not used

