# Synthetic lethality drug repurposing - DepMap and PRISM analysis

This repo supports the manuscript *Challenges and opportunities for drug repurposing in cancers based on synthetic lethality induced by tumor suppressor gene mutations*, under review at npj Systems Biology.

It contains R scripts to build driver-centric dependency networks, run ANOVA- and RF-based tests, integrate drug sensitivity, and generate figures.

## What this repo does

* Builds pan-cancer and tissue-specific mutation→dependency networks from DepMap **CRISPR** data (23Q2).
* Uses **Boruta/Random Forest** to select mutation features linked to gene dependency.
* Adjusts for **MSI** and cancer context (lineage/type) with Type II ANOVA.
* Integrates **PRISM** drug screens (23Q4) via target mapping and correlation tests.
* Exports RDS/CSV summaries and simple plots.

## Repo layout (key scripts)

| File                                    | Purpose                                                                                                                    |
| --------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `functions_RF.R`                | Helper functions: `blank()`, `pull_mut()`, `run_cor_small_RF()`, variants. ANOVA and correlation wrappers.                 |
| `run_feature_select_MUT_DEP_cancer_type_23Q2_MV.R` | Tissue-level loop. Loads precomputed RF hits and runs per-driver analyses within a cancer type. Writes out per-driver RDS. |
| `run_feature_select_MUT_DEP_ALL_23Q2_MV.R` | Pan-cancer loop. Loads precomputed RF hits and runs per-driver analyses. Writes out per-driver RDS. |
| `functions.R`      | Utilities for Cohen’s *d* and basic plotting of KO fitness vs mutation.                                                    |
| `run_cohen_pearson_cancer_type_MV.R`                        | Runs Pearson and Cohen's effect size functions on each genetic interation by cancer type.                                                                 |


## Data inputs (expected)

Set a base path (or edit scripts):

```
DATA_DIR=/global/home/hpc5173/projects/PRISM/data/CCLE
```

**DepMap Public 23Q2** (CSV):

* `CRISPRGeneEffect.csv` 
* `OmicsSomaticMutationsMatrixDamaging.csv`, `OmicsSomaticMutationsMatrixHotspot.csv`
* `OmicsExpressionProteinCodingGenesTPMLogp1.csv`
* `OmicsCNGene.csv`
* `OmicsFusionFiltered.csv`
* `Model.csv`

**DepMap Public 24Q2**:

* `OmicsSignatures.csv` (for MSI, used in ANOVA models)

**PRISM (23Q2/23Q4)**:

* `Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv`
* `Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv`
* `PRISM_d2t_23Q2_citeline_prism.csv` (drug→target mapping; provided)

**TCGA (for mutual exclusivity)**:

* Mutation/CNA matrices per cancer type (see manuscript Methods).

## Software

* R ≥ 4.2 recommended.
* Packages: `tidyverse`, `data.table`, `stringr`, `pbmcapply`/`pbapply`, `Boruta`, `ranger`, `car`, `effsize`, `reshape2`, `ggpubr`.

Install quickly:

```r
install.packages(c(
  "tidyverse","data.table","stringr","pbmcapply","pbapply",
  "Boruta","ranger","car","effsize","reshape2","ggpubr"
))
```

## Quick start

1. **Configure paths.** Edit file paths at the top of scripts or export `DATA_DIR` and keep relative paths.
2. **Preprocess (optional).** Generate any cached blanks (e.g., `blank_model3_*.RDS`) and combined mutation tables you reference.
3. **Run pan-cancer RF feature selection** (example — adjust to your entrypoint):

   ```bash
   Rscript code/MUT_DEP_RF_pan.R /path/to/gene_list.txt
   ```

   Output: `output/RF/MUT_DEPENDENCY_23Q2*/<tissue>/<gene>_getImp*.RDS`
4. **Run tissue-specific ANOVA/summary**:

   ```bash
   Rscript ANOVA_MUT_DEP_byTissue_RF_commented.R BladderCancer
   ```

   Output: `output/RF/COMBINED_MUT_DEP/ANOVA_*/<tissue>/sl_connections_<driver>.RDS`
5. **Cohen’s d and plots (optional):**

   ```r
   source("mut_dependency_utils_commented.R")
   effect_size(target="BRD2", source="ARID1A", type="Damaging")
   p <- plot_fitness_mut_box("BRD2","ARID1A","Damaging", lin=FALSE)
   ```
6. **PRISM integration (example):**
   Use your integration script to join PRISM profiles to genetic hits and compute drug–KO correlations.

### Parallel settings

* Scripts use `pbmcapply` with `mc.cores` set inside code. Pick values that match your HPC quota.

## Methods (short)

* **RF feature selection (Boruta)** on CRISPR effect scores with mutation features (Hotspot + Damaging). Three importance functions: `getImpExtraGini`, `getImpRfZ`, `getImpExtraRaw`. `pValue=0.01`, `maxRuns=500`. A feature is kept only if all three agree.
* **ANOVA (Type II)** with MSI and context covariates. Example models:

  * Pan-cancer: `dep ~ MSI + OncotreeLineage + {dam|hot|dam_b|hot_b}`
  * Tissue_spcific: `dep ~ MSI + {dam|hot|dam_b|hot_b}`
* **Drug integration:** map PRISM drugs to targets; test drug–KO profile correlations and mutation-stratified effects (Cohen’s *d*, t-test/lm).

## Outputs

* `output/RF/MUT_DEPENDENCY_*/*.RDS`: Boruta stats per driver/target.
* `output/RF/COMBINED_MUT_DEP/ANOVA_*/<tissue>/*.RDS`: per-driver association tables.
* `output/COHENS/DEP/*.RDS`: Cohen’s *d* summaries for selected drivers/targets.
* Optional CSV exports for figures/tables (you can add small exporters).

## Data schema notes

* `Model.csv`: uses `ModelID`, `OncotreePrimaryDisease`, `OncotreeLineage`, `GrowthPattern`.
* Mutation matrices: `V1` is model ID; values 0/1/2 for WT/het/hom.
* CRISPR matrices: first column is model ID; others are gene symbols.
* MSI: from `OmicsSignatures.csv`.

## Repro tips

* Set random seeds in any RF wrappers if you change defaults.
* Column **positions** are used in some helpers (e.g., rename indices 7/8/11/12). If your matrix order differs, switch to name-based renames.
* Some code uses `A1BG` as a sentinel column for `drop_na`; confirm it exists in `effect`.
* Impute missing MSI as 0 only if that matches your analysis plan.

## Known limitations (from the study)

* Drug inhibition rarely matches genetic KO profiles at scale.
* Strong context dependence; many hits are cancer type–specific.
* Limited power for rare drivers or small tissue cohorts.



## License

Add a LICENSE file. MIT is common for analysis code.

## Contact

* Maintainer: 0mcv@queensu.ca
* Collaborators: Michael Vermeulen, Andrew Craig, Tomas Babak
