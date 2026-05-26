————————————
GENERAL INFORMATION 
————————————

This readme file was generated on [2026-05-26] by [Emily Meeus] 

Description of Dataset: Dataset used in the paper “ATP13A4 gates extracellular polyamine levels to control excitatory synaptogenesis”. 
This repository contains the R scripts and source data tables used to generate the expression, morphology, and developmental phenotype analyses reported in the paper.

————————
FOLDER STRUCTURE
————————

van-Veen-et-al.-2026/
   
    LICENSE
    README.md

    expression analysis/
        Bakken et al., 2021.R
        Jorstad et al., 2023.R
        Saunders et al., 2018.R
        Schaum et al., 2018.R
        Zhang et al., 2014.R
        Zhang et al., 2016.R

    sholl analysis/
        exp1KO.csv … exp4WT.csv
        key.csv
        script

    Fig. 6 analysis/
        eye_opening.xlsx
        eye_opening_females.xlsx
        eye_opening_males.xlsx
        weight.xlsx
        weight_females.xlsx
        weight_males.xlsx
        script_eyeopening.R
        script_weight.R          

————————
SOFTWARE REQUIREMENTS
————————

- R (version 4.4.3)
- RStudio 

Install once before running:

Data wrangling: install.packages(c("tidyverse", "data.table", "readxl")

Mixed models and stat: install.packages(c("lme4", "lmerTest", "ordinal", "nlme", "car", "multcomp")

Plotting and reshaping: install.packages(c("pheatmap", "reshape")

Rstudio helper: install.packages(c("rstudioapi")

————————
INPUT DATA
————————

Two categories:

1. Included in this repository — the inputs (Excel/CSV) required for the Sholl and Figure 6 analyses are versioned alongside the scripts.
2. Downloaded separately — the expression analyses publicly available datasets. Each script reads files from a working directory the user sets at the top.

### Public datasets used in `expression analysis/`

| Script | Source dataset | Where to obtain |
|---|---|---|
| Bakken et al., 2021.R | Comparative LGN single-nucleus RNA-seq (human, macaque, mouse) | Allen Brain Map — Cell Types Database (LGN, 2021) |
| Jorstad et al., 2023.R | Human MTG single-nucleus RNA-seq | Allen Brain Map (MTG SMART-seq, 2018-06-14 release) |
| Saunders et al., 2018.R | Mouse brain cell atlas (Drop-seq metacells) | DropViz — `metacells.BrainCellAtlas_Saunders_version_2018.04.01.rds` and `annotation...rds` |
| Schaum et al., 2018.R | Tabula Muris (Brain FACS counts + annotations) | figshare / tabula-muris.ds.czbiohub.org |
| Zhang et al., 2014.R | Mouse cortex purified cell-type bulk RNA-seq (FPKM) | brainrnaseq.org |
| Zhang et al., 2016.R | Human cortex purified cell-type bulk RNA-seq | brainrnaseq.org |

Edit the `setwd()` line at the top of each script to point to your local working directory, and place the downloaded files in that working directory.

### Inputs included for `sholl analysis/`

- `exp1KO.csv … exp4WT.csv` — per-experiment Sholl intersection counts. First column = `Radius` (µm); subsequent columns = individual cells/images. Semicolon-separated.
- `key.csv` — maps each input file name to its `Condition` (WT/KO) and replicate number.

### Inputs included for `Fig. 6 analysis/`

- `eye_opening*.xlsx` — columns: `animal`, `genotype` (WT/KO), `day` (postnatal day), `score` (0 = closed, 1 = partial, 2 = fully open).
- `weight*.xlsx` — columns: `animal`, `genotype` (WT/KO), `day` (postnatal day), `weight` (g).
- `*_females` / `*_males` files contain the corresponding sex-stratified subsets.

————————
WORKING EXAMPLE
————————

### Expression analysis

Each `expression analysis/*.R` script is independent.

1. Download the relevant public dataset listed above.
2. Open the script and edit `setwd("...")` at the top to point to the folder containing the dataset.
3. Run the full script. 

### Sholl analysis

1. Open `sholl analysis/script` in RStudio.
2. When prompted by `selectDirectory()`, choose the `sholl analysis/` folder of this repository.
3. Run the script. 

### Fig. 6 analysis

Two scripts, one per phenotype:

- **Eye opening** (`script_eyeopening.R`): cumulative-link mixed model (`ordinal::clmm`). Run on the pooled (`eye_opening.xlsx`) or sex-stratified files.
- **Body weight** (`script_weight.R`): linear mixed model (`lme4::lmer`). Run on the pooled (`weight.xlsx`) or sex-stratified files.

For both: open the script, run line-by-line, and select the `Fig. 6 analysis/` folder when `selectDirectory()` prompts.

————————
LICENCE
————————

This code is released under the **MIT License** — see [`LICENSE`](LICENSE) for the full text.

————————
CITATION
————————

If you use this code, please cite both the accompanying publication (https://www.medrxiv.org/content/10.1101/2025.04.04.25325117v1) and this code repository (https://github.com/emsm3eus-crypto/van-Veen-et-al.-2026).

————————
ACKNOWLEDGEMENTS
————————

This work was funded by the Fonds voor Wetenschappelijk Onderzoek (FWO, Research Foundation Flanders) (G094219N), the Queen Elisabeth Medical Foundation for Neurosciences and Aligning Science Across Parkinson’s (ASAP-000458, ASAP-020607) through the Michael J. Fox Foundation for Parkinson’s Research (MJFF).

