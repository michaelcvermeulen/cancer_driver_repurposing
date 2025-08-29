#!/global/home/hpc5173/miniconda3/envs/main/bin/Rscript

# -----------------------------------------------------------------------------
# Title: PAN-cancer mutation–drug association scan using PRISM logFC profiles
# Author: Michael Vermeulen
# Purpose: For each driver mutation (damaging & hotspot) and each PRISM drug,
#          test associations between mutation status and PRISM viability logFC.
#          Adjust for tissue, growth pattern, and MSI in an ANOVA.
#
# Inputs (update paths as needed):
#   - DepMap CCLE 23Q2: Model.csv, OmicsSomaticMutationsMatrixDamaging.csv,
#                        OmicsSomaticMutationsMatrixHotspot.csv
#   - PRISM 23Q2: Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv
#                 Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv
#                 PRISM_d2t_23Q2_citeline_prism.csv (drug?targets mapping)
#   - Signatures (24Q2): OmicsSignatures.csv (used for MSIScore)
#   - Driver list: drivers_Bailey.csv (column "Gene")
#   - Precomputed pan-cancer ANOVA hits (RDS) used to subset sources: 
#       Pan_important_DEP_ANOVA_PAN_MARCH22_24.RDS (must contain column `source`)
#
# Outputs:
#   - Per-source RDS in /global/home/hpc5173/projects/SL_network/output/PRISM_MUT/PAN2/
#       <SOURCE>.RDS with per-drug rows and columns:
#         pearson, pearson_R, lm_p, t_test, cohens, anova_p, anova_tissue_p,
#         anova_growth_p, anova_msi_p, n, n_mut (and corresponding *_h for hotspot)
#
# Notes:
#   - This script uses Pearson correlation (binary mut vs logFC), t-test, linear
#     model, Cohen's d (effsize::cohen.d), and Type-II ANOVA (car::Anova).
#   - Tissue filtering keeps OncotreePrimaryDisease groups with =10 models.
#   - It analyzes damaging (DAM) features for all sources; hotspot (HOT) is
#     analyzed only if the source gene exists in the HOT matrix.
#   - Parallelization: pbmclapply (forking). Use Linux/macOS. On Windows this
#     will run sequentially. Tune mc.cores as needed.
#   - Required packages beyond those loaded explicitly: effsize, car, reshape2.
# -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(pbmcapply)
library(effsize)  # needed for cohen.d
library(car)      # needed for Anova(type = "II")
library(reshape2) # using reshape2::melt via namespace below

# Pan-cancer ANOVA hits used to subset candidate sources (must contain `source`)
readRDS("/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Pan_important_DEP_ANOVA_PAN_MARCH22_24.RDS") -> o
# If `o$source` includes suffixes, uncomment to split:
# stringr::str_split(string = o$source, pattern = "_", simplify = TRUE)[,1] -> o$mut
# stringr::str_split(string = o$source, pattern = "_", simplify = TRUE)[,2] -> o$type

# Normalize column names (strip trailing " (…)") for consistency
stringr::str_split(string = colnames(effect), pattern = " \\(", simplify = TRUE)[,1] -> colnames(effect)
stringr::str_split(string = colnames(dep),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(dep)
names(dep)[1] <- "V1"  # make model ID explicit if needed later

# ---- PRISM drug matrices and melt to long format -----------------------------
data.table::fread(file = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/PRISM_d2t_23Q2_citeline_prism.csv") -> targets
# Compound metadata (IDs, names, etc.)
data.table::fread(file = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv") -> pr
# PRISM logFC/well-level matrix (wide); will melt to long
data.table::fread(file = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv") -> mat
reshape2::melt(data = mat) -> prism  # columns: V1 (drug ID), variable (ModelID), value (logFC)

# ---- Model annotations, mutation matrices, and MSI signature -----------------
data.table::fread(file = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/Model.csv") -> model
# Normalize tissue labels to match earlier code
model$OncotreePrimaryDisease %<>% gsub(" ", "", .) %>% gsub("-", "", .)

# Mutation matrices (wide, models as rows)
data.table::fread(file = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsSomaticMutationsMatrixDamaging.csv") -> dam
data.table::fread(file = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsSomaticMutationsMatrixHotspot.csv") -> hot
stringr::str_split(string = colnames(dam), pattern = " \\(", simplify = TRUE)[,1] -> colnames(dam)
stringr::str_split(string = colnames(hot), pattern = " \\(", simplify = TRUE)[,1] -> colnames(hot)

# Signatures (used for MSIScore covariate)
data.table::fread(file = "/global/home/hpc5173/projects/PRISM/data/CCLE/24Q2/OmicsSignatures.csv") -> sig

# Driver gene list
data.table::fread(file = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/drivers_Bailey.csv") -> drivers
drivers$Gene %>% unique() -> drivers

# Optionally restrict to drivers present in mutation matrices and pan hits
# (keep overlap with available DAM/HOT cols and `o$source`)
drivers[drivers %in% c(colnames(dam), colnames(hot))] -> gg
gg[gg %in% o$source] -> gg

# ---- Main analysis loop: per source gene, scan all PRISM compounds -----------
lapply(X = gg, FUN = function(x){
  # Skip if cached
  if (file.exists(paste0("/global/home/hpc5173/projects/SL_network/output/PRISM_MUT/PAN2/", x, ".RDS"))) { return(NULL) }
  message(x)
  # Optional: gene-specific skip
  if (x == "BAP1") { return(NULL) }

  # Per-drug parallel loop
  pbmclapply(X = pr$IDs %>% unique(), mc.cores = 16, FUN = function(y){
    message(y)

    # PRISM long rows for this drug ID `y` (value = logFC, variable = ModelID)
    prism %>% dplyr::filter(V1 == y) -> tmp

    # (Optional) drug?targets rows (unused downstream)
    targets %>% dplyr::filter(IDs == y) -> tmp_
    stringr::str_split(string = tmp_$combined_targets, pattern = ",") %>% unique() %>% unlist() -> tr  # not used

    # Keep tissues with =10 models; derive eligible ModelIDs
    model %>%
      dplyr::group_by(OncotreePrimaryDisease) %>%
      dplyr::summarize(count = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(count)) %>%
      dplyr::filter(count >= 10) %>%
      dplyr::pull(OncotreePrimaryDisease) -> types

    model %>% dplyr::filter(OncotreePrimaryDisease %in% types) %>% dplyr::pull(ModelID) -> IDS

    # ---- Damaging (DAM) analysis --------------------------------------------
    # Join DAM mutation for gene x to PRISM logFC for drug y, add covariates
    dam %>% dplyr::select(V1, dplyr::all_of(x)) %>% dplyr::filter(V1 %in% IDS) -> d
    dplyr::left_join(x = d, y = tmp, by = c("V1" = "variable")) %>% tidyr::drop_na(value, V1.y) -> d
    dplyr::left_join(x = d, y = model %>% dplyr::select(ModelID, OncotreePrimaryDisease, GrowthPattern),
                     by = c("V1" = "ModelID")) -> d
    names(d)[3] <- "drug_name"  # from PRISM melted matrix
    names(d)[4] <- "logFC"      # PRISM response
    names(d)[2] <- "mut"        # mutation count/flag for gene x

    # Binarize mutation (>0 ? 1 else 0)
    ifelse(d$mut > 0, 1, 0) -> d$mut

    # Add MSI signature and drop rows missing it
    dplyr::left_join(x = d, y = sig, by = c("V1")) -> d
    d %>% tidyr::drop_na(MSIScore) -> d

    # Evaluate only if at least 4 mutated samples overall
    if (sum(d$mut) > 3) {
      if (length(d[d$mut == 1, ]$mut) <= 2) {
        # Too few mutated samples (<=2): return NA metrics
        data.frame(
          source = x, drug = y,
          pearson = NA, pearson_R = NA, lm_p = NA, t_test = NA, cohens = NA,
          anova_p = NA, anova_tissue_p = NA, anova_growth_p = NA, anova_msi_p = NA,
          n = NA, n_mut = NA
        ) -> output1
      } else {
        # Pearson corr (binary mut vs logFC)
        cor.test(x = d$mut, y = d$logFC, method = "pearson") -> pearson
        # Linear model and t-test
        lm(d, formula = logFC ~ as.factor(mut)) %>% summary() -> linear
        t.test(logFC ~ mut, data = d) -> t_test
        # Cohen's d (requires effsize)
        effsize::cohen.d(f = as.factor(d$mut), d = d$logFC)$estimate -> cohens
        # Type-II ANOVA with covariates (requires car)
        aov(formula = logFC ~ OncotreePrimaryDisease + GrowthPattern + MSIScore + mut, data = d) -> a1
        car::Anova(mod = a1, type = "II", singular.ok = TRUE) -> m1

        data.frame(
          source = x, drug = y,
          pearson = pearson$p.value, pearson_R = pearson$estimate,
          lm_p = linear$coefficients[2, 4], t_test = t_test$p.value,
          cohens = cohens,
          anova_p = m1$`Pr(>F)`[4],
          anova_tissue_p = m1$`Pr(>F)`[1],
          anova_growth_p = m1$`Pr(>F)`[2],
          anova_msi_p = m1$`Pr(>F)`[3],
          n = nrow(d), n_mut = nrow(d[d$mut == 1, ])
        ) -> output1
      }
    } else {
      # Not enough mutated samples total
      data.frame(
        source = x, drug = y,
        pearson = NA, pearson_R = NA, lm_p = NA, t_test = NA, cohens = NA,
        anova_p = NA, anova_tissue_p = NA, anova_growth_p = NA, anova_msi_p = NA,
        n = NA, n_mut = NA
      ) -> output1
    }

    # ---- Hotspot (HOT) analysis (optional) -----------------------------------
    if (x %in% colnames(hot)) {
      hot %>% dplyr::select(V1, dplyr::all_of(x)) %>% dplyr::filter(V1 %in% IDS) -> d
      dplyr::left_join(x = d, y = tmp, by = c("V1" = "variable")) %>% tidyr::drop_na(value, V1.y) -> d
      dplyr::left_join(x = d, y = model %>% dplyr::select(ModelID, OncotreePrimaryDisease, GrowthPattern),
                       by = c("V1" = "ModelID")) -> d
      names(d)[3] <- "drug_name"
      names(d)[4] <- "logFC"
      names(d)[2] <- "mut"
      ifelse(d$mut > 0, 1, 0) -> d$mut

      dplyr::left_join(x = d, y = sig, by = c("V1")) -> d
      d %>% tidyr::drop_na(MSIScore) -> d

      if (length(d[d$mut == 1, ]$mut) <= 2) {
        # NOTE: The next two fields reference `pearson` from the DAM branch.
        # If DAM didn't compute `pearson`, this will error. Kept as-is.
        data.frame(
          source_h = x, drug_h = y,
          pearson_h = pearson$p.value,      # may be undefined
          pearson_R_h = pearson$estimate,   # may be undefined
          lm_p_h = NA, t_test_h = NA, cohens_h = NA,
          anova_p_h = NA, anova_tissue_p_h = NA,
          anova_growth_p_h = NA, anova_msi_p_h = NA,
          n_h = NA, n_mut_h = NA
        ) -> output2
      } else {
        cor.test(x = d$mut, y = d$logFC, method = "pearson") -> pearson
        lm(d, formula = logFC ~ as.factor(mut)) %>% summary() -> linear
        t.test(logFC ~ mut, data = d) -> t_test
        effsize::cohen.d(f = as.factor(d$mut), d = d$logFC)$estimate -> cohens
        aov(formula = logFC ~ OncotreePrimaryDisease + GrowthPattern + MSIScore + mut, data = d) -> a1
        car::Anova(mod = a1, type = "II", singular.ok = TRUE) -> m1

        data.frame(
          source_h = x, drug_h = y,
          pearson_h = pearson$p.value, pearson_R_h = pearson$estimate,
          lm_p_h = linear$coefficients[2, 4], t_test_h = t_test$p.value,
          cohens_h = cohens,
          anova_p_h = m1$`Pr(>F)`[4],
          anova_tissue_p_h = m1$`Pr(>F)`[1],
          anova_growth_p_h = m1$`Pr(>F)`[2],
          anova_msi_p_h = m1$`Pr(>F)`[3],
          n_h = nrow(d), n_mut_h = nrow(d[d$mut == 1, ])
        ) -> output2
      }
    } else {
      # No HOT column for this gene
      data.frame(
        source_h = NA, drug_h = NA,
        pearson_h = NA, pearson_R_h = NA, lm_p_h = NA, t_test_h = NA, cohens_h = NA,
        anova_p_h = NA, anova_growth_p_h = NA, anova_msi_p_h = NA, anova_tissue_p_h = NA,
        n_h = NA, n_mut_h = NA
      ) -> output2
    }

    # Combine DAM and HOT results for this (x, y)
    cbind(output1, output2 %>% dplyr::select(-1, -2)) %>% return()
  }) -> l

  # Bind all drugs for this source gene and attach compound metadata
  data.table::rbindlist(l, fill = TRUE) -> l
  dplyr::left_join(x = l, y = pr, by = c("drug" = "IDs")) %>%
    dplyr::distinct(drug, .keep_all = TRUE) -> l

  # Cache per-source results
  saveRDS(file = paste0("/global/home/hpc5173/projects/SL_network/output/PRISM_MUT/PAN2/", x, ".RDS"), object = l)
  gc()
  return(l)
}) -> ll

# Collate all sources into a single data.table (optional downstream use)
data.table::rbindlist(ll) -> ll

# -----------------------------------------------------------------------------
# Example usage:
#   Rscript PRISM_mut_PAN_commented.R
# Notes:
#   - This script currently reads the pan-hits RDS to subset candidate sources.
#     If you want to drive from a specific gene list, replace `gg` accordingly.
#   - Review the HOT branch note about `pearson` scope before production runs.
# -----------------------------------------------------------------------------
