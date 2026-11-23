# EC-durable-arrest

Network, code and result tables for "Exposure Time, Intervention Sequence and Persistent Signalling Shape Durable
Growth Arrest in a Logical Model of Endometrial Cancer".

Xiaowei Liu, Ning Yu, Xin He, Kehu Yang, Yongxiu Yang

This repository reproduces every number, figure and table of the manuscript from the inputs in `data/`. The literature
evidence for each of the 168 interactions of the network is given in Supplementary Table S4 of the article.

## Quick start

```bash
pip install -r requirements.txt   # or: conda env create -f environment.yml
python run_all.py --figures    # redraw the figures from the result tables in results/tables (under 1 minute)
python run_all.py --new        # recompute scripts 26, 27 and 30-32, then redraw the figures
python run_all.py              # recompute everything (about 1 hour on two CPU cores)
```

Python 3.11 or later. `cmaboss` (MaBoSS) and `mpbn` are needed only for Python scripts 01 to 03; on Windows,
`conda install -c colomoto cmaboss` is the easiest way to install MaBoSS. Each script can also be run on its own.

## Folder structure

```
run_all.py                 runs every Python step in order
bioinformatics/            R scripts of the transcriptomic selection of the four input genes (limma, Cox models)
data/
  EC_network.bnet          the 37-node logical network (rules; Supplementary Table S3)
  signature/               training and testing cohort tables of the candidate genes (Figure 1, Table S1)
  TCGA_UCEC_*.csv          input context, risk group, clinical data, driver mutations, network-gene expression and
                           Hallmark scores of the TCGA tumours
  CPTAC_UCEC_95_*.csv      expression, protein and pathway read-outs of the 95 CPTAC tumours
  GSE236228_*.xlsx, GSE61481_*.xlsx
                           agreement of drug-treated endometrial cancer cell samples with the 32 fixed points
  second_network/          reduced T-LGL leukaemia survival network (six nodes; Saadatpour et al. 2011)
  published_experiments/   Daigh2024_p16 (numerical source data) and Paulsohn2026_withdrawal (digitised curves)
code/
  lib/                     network model, exact Markov-chain engine, figure style,
                           network drawing with horizontal and vertical connectors (orthonet.py, netdraw.py)
  analysis/                analysis scripts; each writes tables to results/tables
  figures/                 one script per figure; writes results/figures (figure1_panels.py redraws panels 1a and 1f)
results/
  tables/                  result tables
  figures/                 one PDF per main-text panel (fig1a.pdf ... fig7i.pdf) and FigureS1-FigureS4 (PDF and 600 dpi PNG)
```

## Analysis scripts

| Script | What it does | Tables | Manuscript |
|---|---|---|---|
| `01_attractor_repertoire.py` | attractors of the 16 input contexts (minimal trap spaces, exact asynchronous dynamics) | A1, A2 | Section 2.2, Fig. S1 |
| `02_context_phenotype_probabilities.py` | probability of reaching each attractor (MaBoSS) | A3 | Fig. 2c, Fig. S1 |
| `03_canalisation_and_trajectories.py` | nodes fixed by the inputs | A6 | Fig. 2b |
| `05_native_switch_exact.py` | exact chain of the seven-node switch: exposure-time responses, schedules, random rates | D1-D5 | Fig. 3b, 4b,c |
| `06_time_window_crossover.py` | fine schedule grid and crossover L* | D4b, D4c | Fig. 4b |
| `07_full_network_pulse_screen.py` | 60 transient interventions on the full network in C0101 | D6 | Fig. 3d |
| `08_context_mechanism_scope.py` | p27, p16 and p21 interventions from every proliferative attractor; simulator validation | D7 | Table S5 |
| `09_controllability_atlas.py` | all single-node interventions from all proliferative attractors | E1 | Fig. 5, Table 1 |
| `10_edge_ablation.py` | deletion of each of the 168 interactions | E2 | Fig. 6a, Table S5 |
| `11_return_edge_combinations.py` | subsets of the four inhibitory inputs to p27 | E2b | Fig. 6a, Table S5 |
| `14_post_withdrawal_kinetics.py` | exact time courses during and after exposure | E4 | Fig. 4a |
| `21_p16_pulse_published_data.py` | published p16 induction and removal experiment (Daigh et al. 2024) | `Daigh2024_p16/` | Fig. 7d-g, Table S6 |
| `22_pancreatic_withdrawal_published_data.py` | published palbociclib withdrawal experiment (Paulsohn et al. 2026) | `Paulsohn2026_withdrawal/` | Fig. 7h,i, Table S6 |
| `23_robustness_checks.py` | transition rates, synchronous updating, random deletion of interactions | R1 | Table S5 |
| `24_driver_lesions.py` | permanent lesions (AKT2 ON, MYC ON, p21 OFF, AKT2 OFF); lesion frequencies in TCGA | F1, F2 | Fig. 6b,c |
| `25_signature_check.py` | repeats the stepwise Cox selection of the four genes in Python | G1 | Fig. 1b, Table S1 |
| `26_transcriptomic_inputs.py` | risk score, survival by risk group, Hallmark scores, risk group by proliferative context | L1-L5 | Fig. 1, Table S1 |
| `27_endometrial_concordance.py` | drug-treated endometrial cancer cells (GSE236228, GSE61481) and CPTAC proteomes | M1-M3 | Fig. 7a-c, Table S6 |
| `30_closed_forms.py` | closed-form durable arrest probabilities in exact rational arithmetic | H1, H2 | Eqs. (1), (2), Appendix A |
| `31_withdrawal_measures.py` | durable arrest probability, L50 and L90, commitment time courses, schedules, lesions, tumour counts | K1-K9 | Figs. 2-6, Table 1 |
| `32_second_network_tlgl.py` | the same analysis on the reduced T-LGL network | J1-J3 | Fig. S4, Table S7 |

`code/lib/exact_engine.py` holds the exact engine for transient interventions. In the code a transient intervention is
called a pulse, and the durable arrest probability is `P_arrest_FP`.

## Bioinformatic part (R)

The four scripts in `bioinformatics/` select the four input genes from the TCGA-UCEC transcriptome. They are not called
by `run_all.py`, because scripts 01 to 03 need the full expression matrix and the univariate Cox table, which are
derived from the public TCGA data. Their outputs, the two cohort tables, are in `data/signature/`.

| Script | What it does |
|---|---|
| `01_differential_expression.R` | tumour versus normal endometrium with limma (absolute log2 fold change above 1, adjusted P below 0.05) |
| `02_univariate_cox_volcano.R` | univariate Cox screen of the immune-related genes (P below 0.05) |
| `03_penalised_cox.R` | stratified 7:3 split (seed 12345) and penalised Cox regression with 10-fold cross-validation |
| `04_stepwise_cox.R` | bidirectional stepwise Cox regression (BIC), risk score, cutoff and performance |

R packages: limma, pheatmap, ggplot2, ggrepel, dplyr, data.table, glmnet, survival, MASS, timeROC.

## Reproducibility

Every random step uses a fixed seed and `run_all.py` fixes `PYTHONHASHSEED`. The exact results (tables D, E, F, H, J, K)
are reproduced exactly. The MaBoSS simulations of scripts 02 and 03 are not exactly reproducible, so basin
weights change in the third decimal place on a rerun; the tables in `results/tables` are those reported.

## Data sources

TCGA-UCEC (Genomic Data Commons; TCGA PanCancer Atlas through cBioPortal), CPTAC-UCEC (cBioPortal), GEO GSE236228
and GSE61481, and the source data of
Daigh et al., Nat. Commun. 2024, 15, 9181. Growth curves of Paulsohn et al., Cell Death Dis. 2026, 17, 129
(doi 10.1038/s41419-025-08362-w) were digitised from the published figure.

## How to cite

Please cite the article and this repository (see `CITATION.cff`; GitHub shows a "Cite this repository" button).

## Licence

Code: MIT (see `LICENSE`). The TCGA-derived tables and the published source data remain under the terms of their sources.
