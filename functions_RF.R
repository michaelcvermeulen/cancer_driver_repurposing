# -----------------------------------------------------------------------------
# Title: Helper functions for mutation–dependency association scans (RF/ANOVA)
# Author: Michael Vermeulen
# Purpose: Provide utility functions used across pipelines to (i) load cached
#          blank ANOVA tables, (ii) pull per-model mutation status for a driver
#          gene, and (iii) compute association statistics between driver
#          mutation features (DAM/HOT; 0/1/2 or binary) and CRISPR gene
#          dependency for a set of target genes.
#
# Functions:
#   - blank(x):              load precomputed placeholder ("blank") ANOVA tables
#                            for cases with insufficient mutated samples.
#   - pull_mut(g, exp, gene_info, table): return per-model mutation status for
#                            driver gene g using a combined mutation table.
#   - run_cor_small_RF(...): adjust for MSI and lineage; includes four models:
#                            DAM, HOT, DAM_binary (dam_b), HOT_binary (hot_b).
#   - run_cor_small_RF_type(...): same as above but models exclude lineage term
#                            (keeps MSI only); useful for subtype-specific scans.
#   - run_cor_small_RF_2(...): correlation-only (Pearson) variant without ANOVA.
#
# Inputs (expected schemas):
#   - effect:   data.frame with CRISPR dependency columns per target gene;
#               first column should be ModelID (named "ModelID").
#   - mut:      per-model mutation table for a single driver; must include
#               columns: DepMap_ID, Hugo_Symbol, and regression features
#               (dam, hot, dam_b, hot_b) after joins/renames below.
#   - driver_status: sample metadata with MSI column; **assumed** to have
#               ModelID in the 1st column and MSI in the 7th (brittle). Consider
#               replacing `select(1,7)` with explicit names for robustness.
#   - genes:    character vector of target gene symbols present in `effect`.
#
# Notes:
#   - Parallelization uses pbmclapply (forking): Linux/macOS only by default.
#   - `car::Anova(type = "II")` is used; ensure `library(car)` is attached in
#     the calling script. The same for `pbmcapply` and `data.table`.
#   - Column renames by numeric position (7?dam, 8?hot, 11?dam_b, 12?hot_b)
#     assume a stable column order after joins; this is fragile. If schema
#     changes, prefer renaming by column names instead of index.
#   - `tidyr::drop_na(tt, A1BG)` uses A1BG as a sentinel column to drop rows not
#     present in the dependency matrix. Verify that A1BG exists in `effect`.
# -----------------------------------------------------------------------------

blank <- function(x){
  # Load a cached placeholder ANOVA table for the given feature type `x`.
  # Expected `x` values include: "dam", "hot", "dam_b", "hot_b".
  # The RDS should contain a data.frame/matrix with rownames including
  #   c("MSI", "OncotreeLineage", x) so that indexing like m1["MSI", ] works.
  readRDS(paste0("/global/home/hpc5173/projects/SL_network/data/blank_model3_", x, ".RDS")) -> xx
  return(xx)
}



pull_mut <- function(g, exp, gene_info, table){
  # Return per-model mutation status rows for driver gene `g`.
  # - `table` is the combined mutation table (long or pre-shaped) with columns
  #    including `Hugo_Symbol` and `DepMap_ID`.
  # - `gene_info` is used here to annotate chromosome (first column assumed to be
  #    chromosome; brittle). Consider selecting by explicit column name.
  
  table %>% dplyr::filter(Hugo_Symbol == g) -> df
  if (nrow(df) == 0) { return(NULL) }
  
  gene_info %>% dplyr::filter(symbol == g) %>% .[1,1] %>% as.character() -> chr
  df$chr <- chr

  rownames(df) <- df$DepMap_ID
  return(df)
}


run_cor_small_RF <- function(mut, driver_status, effect, genes, all = F){
  # Association pipeline with lineage and MSI covariates + Type-II ANOVA.
  # For each target in `genes`, fit:
  #   dep ~ MSI + OncotreeLineage + {dam|hot|dam_b|hot_b}
  # and compute Pearson correlation between the mutation feature and dependency.

  # Load placeholder ANOVA tables used when mutated sample counts are too small.
  blank("dam")    -> blank_dam
  blank("hot")    -> blank_hot
  blank("dam_b")  -> blank_dam_b
  blank("hot_b")  -> blank_hot_b

  # Validate/limit target list to those present in `effect` (exclude ModelID col)
  if (isTRUE(all)) {
    genes <- colnames(effect) %>% unique() %>% .[-1]
  } else {
    genes[genes %in% (colnames(effect) %>% unique() %>% .[-1])] -> genes
  }
  if (length(genes) == 0) { return(NULL) }

  # Keep the driver gene symbol for output
  mut$Hugo_Symbol %>% unique() %>% as.character() -> mut_gene

  # Join mutation?dependency, then attach MSI from driver_status
  dplyr::left_join(y = effect, x = mut, by = c("DepMap_ID" = "ModelID")) -> tt
  dplyr::left_join(y = driver_status %>% dplyr::select(1,7),  # assumes (ModelID, MSI)
                   x = tt, by = c("DepMap_ID" = "ModelID")) -> tt

  # Require presence of a sentinel gene column (A1BG) to keep rows
  tidyr::drop_na(tt, A1BG) -> tt

  pbmcapply::pbmclapply(X = genes, mc.cores = 8, FUN = function(y){
    # Copy working frame and rename columns by position (fragile)
    tt -> df
    colnames(df)[which(colnames(df) == y)] <- "dep"
    colnames(df)[7]  <- "dam"
    colnames(df)[8]  <- "hot"
    colnames(df)[11] <- "dam_b"
    colnames(df)[12] <- "hot_b"

    # Optionally exclude very small lineages (hard-coded list)
    df %>% dplyr::filter(!(OncotreeLineage %in% c("Ampulla of Vater","Testis","Vulva/Vagina","Adrenal Gland"))) -> df

    # Intended NA fill for MSI (see note above)
    if (!any(is.na(df$MSI))) { df[which(is.na(df$MSI)), ]$MSI <- 0 }

    # DAM model (require >2 mutated)
    if ((sum(df[,7] > 0)) > 2) {
      aov(formula = dep ~ MSI + OncotreeLineage + dam, data = df) -> a1
      l1 <- cor(x = df$dam, y = df$dep, method = "pearson")
      car::Anova(mod = a1, type = "II", singular.ok = TRUE) -> m1
    } else { m1 <- blank_dam; l1 <- NA }

    # HOT model
    if ((sum(df[,8] > 0)) > 2) {
      aov(formula = dep ~ MSI + OncotreeLineage + hot, data = df) -> a2
      l2 <- cor(x = df$hot, y = df$dep, method = "pearson")
      car::Anova(mod = a2, type = "II", singular.ok = TRUE) -> m2
    } else { m2 <- blank_hot; l2 <- NA }

    # DAM_burden model
    if ((sum(df[,11] > 0)) > 2) {
      aov(formula = dep ~ MSI + OncotreeLineage + dam_b, data = df) -> a3
      l3 <- cor(x = df$dam_b, y = df$dep, method = "pearson")
      car::Anova(mod = a3, type = "II", singular.ok = TRUE) -> m3
    } else { m3 <- blank_dam_b; l3 <- NA }

    # HOT_burden model
    if ((sum(df[,12] > 0)) > 2) {
      aov(formula = dep ~ MSI + OncotreeLineage + hot_b, data = df) -> a4
      l4 <- cor(x = df$hot_b, y = df$dep, method = "pearson")
      car::Anova(mod = a4, type = "II", singular.ok = TRUE) -> m4
    } else { m4 <- blank_hot_b; l4 <- NA }

    # Collate statistics into one row per target gene
    data.table::data.table(
      mut = mut_gene,
      target = y,

      MSI_p_dam = m1["MSI", ]$`Pr(>F)`,
      MSI_F_dam = m1["MSI", ]$`F value`,
      OncotreeLineage_p_dam = m1["OncotreeLineage", ]$`Pr(>F)`,
      OncotreeLineage_F_dam = m1["OncotreeLineage", ]$`F value`,
      dam_p  = m1["dam", ]$`Pr(>F)`,
      dam_F  = m1["dam", ]$`F value`,
      dam_R  = l1,

      MSI_p_hot = m2["MSI", ]$`Pr(>F)`,
      MSI_F_hot = m2["MSI", ]$`F value`,
      lineage_p_hot = m2["OncotreeLineage", ]$`Pr(>F)`,
      lineage_F_hot = m2["OncotreeLineage", ]$`F value`,
      hot_p  = m2["hot", ]$`Pr(>F)`,
      hot_F  = m2["hot", ]$`F value`,
      hot_R  = l2,

      MSI_p_dam_b = m3["MSI", ]$`Pr(>F)`,
      MSI_F_dam_b = m3["MSI", ]$`F value`,
      lineage_p_dam_b = m3["OncotreeLineage", ]$`Pr(>F)`,
      lineage_F_dam_b = m3["OncotreeLineage", ]$`F value`,
      dam_b_p  = m3["dam_b", ]$`Pr(>F)`,
      dam_b_F  = m3["dam_b", ]$`F value`,
      dam_b_R  = l3,

      MSI_p_hot_b = m4["MSI", ]$`Pr(>F)`,
      MSI_F_hot_b = m4["MSI", ]$`F value`,
      lineage_p_hot_b = m4["OncotreeLineage", ]$`Pr(>F)`,
      lineage_F_hot_b = m4["OncotreeLineage", ]$`F value`,
      hot_b_p  = m4["hot_b", ]$`Pr(>F)`,
      hot_b_F  = m4["hot_b", ]$`F value`,
      hot_b_R  = l4,

      dam_n  = sum(df$dam > 0),  dam_2n = sum(df$dam == 2), dam_1n = sum(df$dam == 1),
      dam_0n = sum(df$dam == 0), hot_2n = sum(df$hot == 2), hot_1n = sum(df$hot == 1),
      hot_n  = sum(df$hot > 0),  hot_0n = sum(df$hot == 0),
      n = nrow(df)
    ) -> dat

    return(dat)
  }) -> l

  data.table::rbindlist(l) -> l
  return(l)
}


run_cor_small_RF_type <- function(mut, driver_status, effect, genes, all = F){
  # Variant without lineage term in the ANOVA; adjusts only for MSI.

  blank("dam")    -> blank_dam
  blank("hot")    -> blank_hot
  blank("dam_b")  -> blank_dam_b
  blank("hot_b")  -> blank_hot_b

  if (isTRUE(all)) {
    genes <- colnames(effect) %>% unique() %>% .[-1]
  } else {
    genes[genes %in% (colnames(effect) %>% unique() %>% .[-1])] -> genes
  }
  if (length(genes) == 0) { return(NULL) }

  mut$Hugo_Symbol %>% unique() %>% as.character() -> mut_gene

  dplyr::left_join(y = effect, x = mut, by = c("DepMap_ID" = "ModelID")) -> tt
  dplyr::left_join(y = driver_status %>% dplyr::select(1,7),
                   x = tt, by = c("DepMap_ID" = "ModelID")) -> tt

  tidyr::drop_na(tt, A1BG) -> tt

  pbmcapply::pbmclapply(X = genes, mc.cores = 8, FUN = function(y){
    tt -> df
    colnames(df)[which(colnames(df) == y)] <- "dep"
    colnames(df)[7]  <- "dam"
    colnames(df)[8]  <- "hot"
    colnames(df)[11] <- "dam_b"
    colnames(df)[12] <- "hot_b"

    # Remove rare lineages (kept from original; though lineage not modeled here)
    df %>% dplyr::filter(!(OncotreeLineage %in% c("Ampulla of Vater","Testis","Vulva/Vagina","Adrenal Gland"))) -> df

    # Intended MSI NA fill (see note in header)
    if (!any(is.na(df$MSI))) { df[which(is.na(df$MSI)), ]$MSI <- 0 }

    if ((sum(df[,7] > 0)) > 2) {
      aov(formula = dep ~ MSI + dam, data = df) -> a1
      l1 <- cor(x = df$dam, y = df$dep, method = "pearson")
      car::Anova(mod = a1, type = "II", singular.ok = TRUE) -> m1
    } else { m1 <- blank_dam; l1 <- NA }

    if ((sum(df[,8] > 0)) > 2) {
      aov(formula = dep ~ MSI + hot, data = df) -> a2
      l2 <- cor(x = df$hot, y = df$dep, method = "pearson")
      car::Anova(mod = a2, type = "II", singular.ok = TRUE) -> m2
    } else { m2 <- blank_hot; l2 <- NA }

    if ((sum(df[,11] > 0)) > 2) {
      aov(formula = dep ~ MSI + dam_b, data = df) -> a3
      l3 <- cor(x = df$dam_b, y = df$dep, method = "pearson")
      car::Anova(mod = a3, type = "II", singular.ok = TRUE) -> m3
    } else { m3 <- blank_dam_b; l3 <- NA }

    if ((sum(df[,12] > 0)) > 2) {
      aov(formula = dep ~ MSI + hot_b, data = df) -> a4
      l4 <- cor(x = df$hot_b, y = df$dep, method = "pearson")
      car::Anova(mod = a4, type = "II", singular.ok = TRUE) -> m4
    } else { m4 <- blank_hot_b; l4 <- NA }

    data.table::data.table(
      mut = mut_gene,
      target = y,

      MSI_p_dam = m1["MSI", ]$`Pr(>F)`,
      MSI_F_dam = m1["MSI", ]$`F value`,
      dam_p  = m1["dam", ]$`Pr(>F)`,
      dam_F  = m1["dam", ]$`F value`,
      dam_R  = l1,

      MSI_p_hot = m2["MSI", ]$`Pr(>F)`,
      MSI_F_hot = m2["MSI", ]$`F value`,
      hot_p  = m2["hot", ]$`Pr(>F)`,
      hot_F  = m2["hot", ]$`F value`,
      hot_R  = l2,

      MSI_p_dam_b = m3["MSI", ]$`Pr(>F)`,
      MSI_F_dam_b = m3["MSI", ]$`F value`,
      dam_b_p  = m3["dam_b", ]$`Pr(>F)`,
      dam_b_F  = m3["dam_b", ]$`F value`,
      dam_b_R  = l3,

      MSI_p_hot_b = m4["MSI", ]$`Pr(>F)`,
      MSI_F_hot_b = m4["MSI", ]$`F value`,
      hot_b_p  = m4["hot_b", ]$`Pr(>F)`,
      hot_b_F  = m4["hot_b", ]$`F value`,
      hot_b_R  = l4,

      dam_n  = sum(df$dam > 0),  dam_2n = sum(df$dam == 2), dam_1n = sum(df$dam == 1),
      dam_0n = sum(df$dam == 0), hot_2n = sum(df$hot == 2), hot_1n = sum(df$hot == 1),
      hot_n  = sum(df$hot > 0),  hot_0n = sum(df$hot == 0),
      n = nrow(df)
    ) -> dat

    return(dat)
  }) -> l

  data.table::rbindlist(l) -> l
  return(l)
}


