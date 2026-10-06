# MALAT1/miR-204-5p/SIRT1 Boolean Model of Pyroptosis in HCC

Code and model files accompanying the manuscript:

> **Dynamic Modeling Reveals lncRNA MALAT1-mediated Suppression of GSDMD-dependent Pyroptosis through the miR-204-5p/SIRT1 Axis in HCC**
> Epitácio Farias, Laís de Carvalho Gonçalves, Bruna H. C. D. Barbosa, Maria Eduarda S. Pinheiro, Gabriel Rolim, Otávio P. C. de Sousa, Shantanu Gupta\*
> Bioinformatics Multidisciplinary Environment (BioME), Metropolis Digital Institute, UFRN, Natal, Brazil

## Overview

Hepatocellular carcinoma (HCC) is marked by resistance to cell death. This project uses a literature-curated **Boolean network model** to ask whether the lncRNA **MALAT1 / miR-204-5p / SIRT1** axis, already known to drive HCC progression and therapy resistance, also controls **GSDMD-dependent pyroptosis**, an inflammatory and immunogenic form of cell death.

The model (31 nodes, 69 interactions) integrates DNA damage response signaling (input node `DDR`), p53 activity states, the non-coding RNA layer (MALAT1, miR-204-5p), SIRT1, and downstream regulators of five phenotypes:

`PROLIFERATION` · `RESISTANCE` · `CELL_CYCLE_ARREST` · `APOPTOSIS` · `PYROPTOSIS`

Main findings:

- **Wild type:** five stable states. Without DNA damage the network settles in proliferation + resistance; with DNA damage it reaches resistance, cell-cycle arrest, apoptosis or pyroptosis.
- **Perturbations:** MALAT1 overexpression, miR-204-5p knockout and SIRT1 overexpression all lead to resistance; miR-204-5p overexpression and SIRT1 knockout unlock cell-death states, including pyroptosis.
- **Feedback loops:** two double-negative (positive) loops, `miR-204-5p ⊣ SIRT1 ⊣ p53 → miR-204-5p` and `p53 ⊣ MALAT1 ⊣ p53`, act as switches between resistance and cell death.
- **Patient data:** expression of the network nodes was compared between tumor and adjacent normal tissue in **TCGA-LIHC** as an independent, population-level layer of support.

All network-level predictions are computational and intended as testable hypotheses.

## Repository contents

| File | Language | Description |
|---|---|---|
| `MODEL_mir204_pyroptosis_hcc.zginml` | GINsim | The Boolean model (31 nodes, 69 interactions) with its logical rules. Opens directly in GINsim. |
| `Figure2_WT_Stables_States_FIXED.ipynb` | Python | Builds Figure 2: node-activity matrix of the wild-type stable states and their probability of convergence, parsed from the GINsim Monte Carlo output (`.csv`). |
| `code_3_figure.R` | R | Imports the GINsim model into BoolNet, enumerates the exact fixed points with DDR = ON for the wild type and each perturbation (MALAT1 E1, miR-204-5p E1/KO, SIRT1 E1/KO), and builds Figure 3 (model predictions vs. literature). |
| `tcga_lihc_network_analysis.Rmd` | R | TCGA-LIHC analysis: data download, clinical table integration, normalization, tumor vs. normal comparison of the network nodes (Figure 5), and stratification by radiotherapy/chemotherapy. |

## Requirements

**Boolean model**

- [GINsim](http://ginsim.org/) 3.0.0b (logical model construction, stable-state search, Monte Carlo simulations, perturbations and circuit analysis)

**Python** (figure notebook)

- Python ≥ 3.8
- `matplotlib`, `pandas`, `pillow`, `jupyter`

```bash
pip install matplotlib pandas pillow jupyter
```

**R**

- R ≥ 4.2
- Fixed-point analysis (`code_3_figure.R`), CRAN: `BoolNet`, `ggplot2`, `xml2`
- TCGA analysis (`tcga_lihc_network_analysis.Rmd`):
  - Bioconductor: `TCGAbiolinks`, `SummarizedExperiment`, `rtracklayer`, `DESeq2`
  - CRAN: `tidyverse`, `ggsignif`, `rmarkdown`
  - GENCODE v50 annotation (`gencode.v50.annotation.gff3`), available from [gencodegenes.org](https://www.gencodegenes.org/human/)

```r
install.packages(c("BiocManager", "BoolNet", "xml2", "tidyverse", "ggsignif", "rmarkdown"))
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "rtracklayer", "DESeq2"))
```

## Usage

### 1. Boolean model simulations

Open `MODEL_mir204_pyroptosis_hcc.zginml` in GINsim and run the stable-state and Monte Carlo analyses. Phenotype probabilities are estimated from asynchronous Monte Carlo simulations with random initial states; perturbations are done by fixing nodes to 0 (knockout, KO) or 1 (ectopic expression, E1). Export the Monte Carlo results as a `.csv` file.

### 2. Stable-state figure (Figure 2)

Open `Figure2_WT_Stables_States_FIXED.ipynb`, set `CSV_PATH` in the **CONFIGURAÇÃO** cell to the GINsim Monte Carlo `.csv`, and run all cells. The panel is saved as `network_panel.png`.

### 3. Exact fixed points under perturbation (Figure 3)

Run the script passing the model file as the first argument:

```bash
Rscript code_3_figure.R MODEL_mir204_pyroptosis_hcc.zginml
```

The script reads the logical rules from the `.zginml`, converts them to BoolNet, fixes `DDR = 1`, and enumerates the exact fixed points (SAT-based search, no Monte Carlo) for the wild type and each perturbation. Results go to `GINsim_DDR_ON_exact_results/`:

- extracted rules and the BoolNet network (`01`–`03`)
- fixed points, condition summary and phenotype reachability tables (`04`–`06`)
- stable-state and reachability heatmaps (`07`–`10`)
- `Figure_3_model_predictions_vs_literature.png/.pdf`
- methodological report and `sessionInfo()` (`11`–`12`)

Fixed points indicate which phenotypes are reachable under each condition. They are not probabilities or cell frequencies.

### 4. TCGA-LIHC expression analysis (Figure 5)

1. Run the `download-data` chunk of `tcga_lihc_network_analysis.Rmd` **manually once** (it is set to `eval=FALSE` because it downloads large files from the GDC). Adjust the path to the GENCODE GFF3 file before running. This creates `counts.RData`, `mirnaseq_counts.RData` and `clin.RData`.
2. Knit the document (`rmarkdown::render("tcga_lihc_network_analysis.Rmd")`) or run the remaining chunks in order.

Main outputs:

- `master_clinical_lihc.csv` — one row per patient, integrating the GDC clinical table and the BCR Biotab supplementary tables
- `wilcoxon_*.csv` — Wilcoxon rank-sum tests with Benjamini-Hochberg correction
- `*_tumor_vs_normal.png` — boxplots of normalized expression (DESeq2 `vst`) per network node

## Methods summary

| Step | Approach |
|---|---|
| Network curation | Experimentally validated interactions from PubMed, checked against BioGRID and TargetScan |
| Dynamics | Asynchronous Boolean updating in GINsim 3.0; attractors interpreted as phenotypes |
| Phenotype probabilities | Monte Carlo simulations from random initial states |
| Perturbations | Node fixation (KO = 0, E1 = 1) and edge deletion within feedback loops |
| Exact fixed points | BoolNet SAT-restricted enumeration with DDR = 1, from the rules exported by GINsim |
| Patient data | TCGA-LIHC RNA-seq and miRNA-seq via TCGAbiolinks; DESeq2 `vst`; Wilcoxon test + BH correction |

## Citation

If you use this model or code, please cite:

```
Farias E, Gonçalves LC, Barbosa BHCD, Pinheiro MES, Rolim G, Sousa OPC, Gupta S.
Dynamic Modeling Reveals lncRNA MALAT1-mediated Suppression of GSDMD-dependent
Pyroptosis through the miR-204-5p/SIRT1 Axis in HCC. (manuscript in preparation)
```

## Acknowledgements

This work started as a research project in the graduate course BIF0050 – *Tópicos Avançados em Bioinformática IV* ("Hands-On Computational Systems Biology: Modeling Cell Death Pathways in Cancer"), BioME/IMD/UFRN. Computational resources were provided by NPAD/UFRN and BioME/UFRN.

## Contact

- Shantanu Gupta (corresponding author) — shantanu.gupta@imd.ufrn.br
- Epitácio Farias — [@epfarias](https://github.com/epfarias)
