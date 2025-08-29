# -----------------------------------------------------------------------------
# Title: Utilities for mutation–CRISPR fitness analyses (DepMap/PRISM 23Q2)
# Author: Michael Vermeulen
# Purpose:
#   - Compute Cohen's d for KO fitness differences between mutated vs WT lines.
#   - Plot KO fitness distributions by mutation status, optionally faceted by
#     lineage.
#
# Inputs (loaded below; update paths as needed):
#   - Damaging mutations: OmicsSomaticMutationsMatrixDamaging.csv (23Q2)
#   - Hotspot   mutations: OmicsSomaticMutationsMatrixHotspot.csv (23Q2)
#   - CRISPR fitness ("effect"): CRISPRGeneEffect.csv (23Q2)
#
# Additional objects expected by functions (NOT loaded here):
#   - `model`: DepMap Model.csv with at least columns [ModelID, OncotreeLineage]
#
# Packages:
#   - effsize (for Cohen's d), tidyverse, data.table, Boruta (not used here),
#     pbapply (not used here), stringr.
#   - For plotting: ggpubr (ggboxplot, ggpar, ggscatter, stat_cor). Attach it in
#     the calling script: `library(ggpubr)`.
#
# Notes:
#   - Column name normalization strips trailing " (.." from input matrices to
#     keep base gene symbols.
#   - `effect` here uses CRISPRGeneEffect (not CRISPRGeneDependency); confirm
#     this is the intended response for your analysis.
#   - Captions in plots reference "DepMap 22Q4" in the original code but data
#     paths are 23Q2. Caption text is kept; change if needed.
# -----------------------------------------------------------------------------

library(tidyverse)
library(Boruta)   
library(data.table)
library(stringr)
library(pbapply)  
library(effsize)

# ---- Load core matrices (23Q2) ----------------------------------------------
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsSomaticMutationsMatrixDamaging.csv", nThread = 8) -> dam
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/OmicsSomaticMutationsMatrixHotspot.csv",  nThread = 8) -> hot
fread(input = "/global/home/hpc5173/projects/PRISM/data/CCLE/23Q2/CRISPRGeneEffect.csv",                     nThread = 8) -> effect

# Normalize column names to base symbols (drop trailing " (…)")
stringr::str_split(string = colnames(dam),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(dam)
stringr::str_split(string = colnames(hot),    pattern = " \\(", simplify = TRUE)[,1] -> colnames(hot)
stringr::str_split(string = colnames(effect), pattern = " \\(", simplify = TRUE)[,1] -> colnames(effect)


# -----------------------------------------------------------------------------
# Function: effect_size
# Description:
#   Compute Cohen's d for KO fitness (from `effect`) between mutated vs WT cell
#   lines for a given source mutation (type = "Hotspot" or "Damaging") and a
#   CRISPR target gene.
# Returns:
#   One-row data.frame with effect size, group sizes, means, and mean difference.
# Edge cases:
#   - If `source` is not present in the chosen mutation matrix, or <3 mutated
#     cell lines exist, returns a row of NAs (except identifying fields).
# -----------------------------------------------------------------------------
effect_size <- function(target, source, type){
  
  # Pre-allocate an empty result row to return on early exits
  data.frame( target = target,
              source = source,
              type = type,
              estimate_cohens_d = NA,
              sd = NA,
              mag = NA,
              n_mut = NA, 
              n_wt = NA,
              mean_mut = NA,
              mean_wt = NA, diff = NA) -> empty
  
  # Select mutation table by requested type
  if(type == "Hotspot"){ mut <- hot }
  if(type == "Damaging"){ mut <- dam }
  
  # If source gene not present, return empty row
  if(!( source %in% colnames(mut))){ return(empty) }
  
  # Build mutation status table: DepMap_ID + mutation count/flag for `source`
  mut %>% dplyr::select(V1, dplyr::all_of(source)) -> tmp
  names(tmp) <- c("DepMap_ID","MUT")
  
  # Partition models into mutated (1/2) vs wild-type (0)
  tmp %>% dplyr::filter(MUT %in% c(1,2)) %>% .[,1] %>% unlist() %>% as.vector() -> mutated
  tmp %>% dplyr::filter(MUT %in% c(0))   %>% .[,1] %>% unlist() %>% as.vector() -> wt
  
  # Require at least 3 mutated lines for a minimally stable estimate
  if(length(mutated) < 3){ return(empty) }
  
  # Collect fitness values for target gene in mutated vs WT groups
  effect %>% dplyr::filter(ModelID %in% mutated) %>%
    dplyr::select(ModelID, dplyr::all_of(target)) %>% .[,2] %>% unlist() %>% as.vector() -> a
  effect %>% dplyr::filter(ModelID %in% wt) %>%
    dplyr::select(ModelID, dplyr::all_of(target)) %>% .[,2] %>% unlist() %>% as.vector() -> b
  
  # Cohen's d with pooled SD
  cohen.d(a, b, pooled = TRUE) -> cohen
  
  # Assemble tidy one-row output
  data.frame(target = target,
             source = source,
             type = type,
             estimate_cohens_d = cohen$estimate,
             sd = cohen$sd,
             mag = cohen$magnitude,
             n_mut = length(a), 
             n_wt = length(b),
             mean_mut = mean(a),
             mean_wt = mean(b)) -> df
  
  df %>% dplyr::mutate(diff = mean_mut - mean_wt ) %>% return()
}


# -----------------------------------------------------------------------------
# Function: plot_fitness_mut_box
# Description:
#   Produce boxplots of KO fitness ("eff") by mutation count (0/1/2) for a given
#   source?target pair. Optionally facet by OncotreeLineage when `lin = TRUE`.
# Requirements:
#   - `model` data.frame in the calling environment with columns [ModelID,
#      OncotreeLineage].
#   - `ggpubr` attached for ggboxplot, ggpar, ggscatter, stat_cor.
# Arguments:
#   - target_gene: CRISPR KO gene symbol in `effect` matrix
#   - source_gene: mutation gene symbol in `dam`/`hot`
#   - t: "Damaging" or "Hotspot" (selects the mutation matrix)
#   - lin: if TRUE, facet by OncotreeLineage
# Returns:
#   A list of two ggplot objects: [boxplot, scatter + stat_cor]
# -----------------------------------------------------------------------------
plot_fitness_mut_box <- function(target_gene, source_gene, t, lin = FALSE){
  
  # Choose mutation matrix
  t %>% as.character() -> t
  if(t == "Damaging"){ mat <- dam }
  if(t == "Hotspot"){ mat <- hot }
  
  # Pull CRISPR fitness for target gene; rename to `eff`
  effect %>% dplyr::select(V1, dplyr::all_of(target_gene)) -> dd
  names(dd)[2] <- "eff"
  
  # Join mutation matrix to fitness by model ID, then attach lineage from `model`
  dplyr::left_join(y = mat, x = dd, by = 'V1') -> df4
  dplyr::left_join(y = df4, x = model %>% dplyr::select(1,24), by = c("ModelID" = "V1")) -> df4
  
  # Keep complete cases and then select just the relevant columns
  df4[complete.cases(df4)] -> df4
  df4 %>% dplyr::select(1:3, dplyr::all_of(source_gene)) -> dff
  
  # Boxplot: by mutation count (0/1/2); optionally facet by lineage
  if(isFALSE(lin)){
    ggboxplot(
      data = dff,
      x = dplyr::all_of(source_gene), y = "eff",
      add = "jitter", add.params = list(shape = 21, fill = NA),
      fill = dplyr::all_of(source_gene), palette = "aaas"
    ) -> p
    ggpar(
      p,
      xlab = paste0("Source gene: ", source_gene, " MUT"),
      ylab = paste0("Target gene: ", target_gene, " KO fitness"),
      legend = "none",
      caption = paste(tools::toTitleCase(t), "mutations (DepMap 22Q4)")  # caption text retained from original
    ) -> p
  } else {
    ggboxplot(
      data = dff,
      x = dplyr::all_of(source_gene), y = "eff",
      add = "jitter", facet.by = 'OncotreeLineage',
      add.params = list(shape = 21, fill = NA),
      fill = dplyr::all_of(source_gene), palette = "aaas"
    ) -> p
    ggpar(
      p,
      xlab = paste0("Source gene: ", source_gene, " MUT"),
      ylab = paste0("Target gene: ", target_gene, " KO fitness"),
      legend = "none",
      caption = paste(tools::toTitleCase(t), "mutations (DepMap 22Q4)")
    ) -> p
  }
  
  # Prepare scatter with regression line and correlation annotation
  names(dff) <- c("Model","Onco","eff","mut")
  ifelse(dff$mut > 0, 1, 0) -> dff$mut_binary
  
  # Quick diagnostics (printed to console)
  cor.test(x = dff$eff, y = dff$mut) -> co
  lm(data = dff, formula = mut ~ eff ) %>% summary() -> s
  print(paste0("t-value: ", round(s$coefficients[2,3], 3)))
  print(paste0("Correlation coef: ", round(co$estimate, 3)))
  
  # Scatter of mutation count vs fitness (with fitted line and r stat)
  dff$mut %>% as.factor() -> dff$mut_fac
  ggscatter(data = dff, x = "mut", y = "eff", add = "reg.line", shape = 21, fill = NA) -> p2
  ggpar(
    p2,
    xlab = paste0("Source gene: ", source_gene, " MUT"),
    ylab = paste0("Target gene: ", target_gene, " KO fitness"),
    legend = "none",
    caption = paste(tools::toTitleCase(t), "mutations (DepMap 22Q4)")
  ) -> p2
  
  return(list(p, p2 + stat_cor()))
}
