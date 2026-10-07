# 🧬 Methylation Aging: Multi-Model Epigenetic Clock Analysis

This repository implements a comparative framework evaluating multiple DNA methylation foundation models and sequence-based models for predicting unsupervised biological age across independent rhesus macaque (*Macaca mulatta*) tissues.

---

## 📂 Project Structure

Tracked source. The `*.csv` / `*.RData` matrices are gitignored and are **not**
in this repository — see [Data](#data).

```
├── prepare_data.py                   # extract liver subset, fetch Ensembl sequence windows
├── cpg_liftover.py                   # CpG coordinate liftover to Mmul_10
├── extract_subset.R                  # pull the top variable sites from the raw RData
├── train_models_multimodel.ipynb     # primary analysis & evaluation notebook
├── train_models_multimodel.py        # script form of the above
├── train_models_cpg_density.py       # CpG-density variant
├── cayo_analysis.ipynb               # Cayo blood cohort
├── cpgpt_analysis.ipynb              # CpGPT foundation-model run (ONPRC liver)
├── selected_window_sequences.json    # DNA sequence windows, blood CpG sites
└── liver_selected_sequences.json     # DNA sequence windows, liver CpG sites
```

---

## 🧬 Project Overview
This project compares traditional regularized regression against representations derived from DNA language models (LLMs) and methylation transformer models. 

Models are trained on 90% of the dataset and evaluated on a **10% holdout test set** to predict an unsupervised biological age proxy calculated directly from global methylation statistics.

---

## 📊 Datasets

### 1. Blood Dataset (Cayo Santiago Macaques)
- **Samples**: 572 longitudinal blood RRBS samples from 216 unique individuals.
- **Features**: 216 CpG window methylation percentages (`selected_percent_methylation.csv`).
- **Splitting**: Evaluation uses a **Grouped Split** by `animal_id` to ensure that no samples from the same animal leak between the training and test sets.

### 2. Liver Dataset (ONPRC Cohort)
- **Samples**: 94 liver RRBS samples.
- **Features**: The top 200 most variable CpG sites selected from a total of 235,810 sites in the raw matrix (`ONPRCliverSitePercMeth5x1090n94.RData`).
- **Sequences**: Genomic coordinates lifted over and mapped to extract `+/- 100bp` context window sequences via the Ensembl REST API for *Macaca mulatta* (assembly `Mmul_10`).

### Upstream pipeline

The methylation matrices consumed here are produced by
[`RheMacAge`](https://github.com/elisabethgoldman/RheMacAge) — the RRBS
preprocessing and clock-construction pipeline (Bismark genome prep → trim →
paired-end align → methylation extraction → 1 kb windowing → coverage filtering
→ imputation → elastic-net clock).

Two things that matter for anything built downstream of this repo:

- The `chrN_start_end` column keys used throughout (`chr4_56044000_56045000`)
  are produced by that pipeline's **`mapWindows.sh`**, which bins CpG-site
  coverage files into 1 kb windows. Any project joining onto these features
  should adopt the same binning rather than re-deriving it.
- Windows are retained only where covered in ≥90% of samples, then imputed and
  filtered at <5× median coverage and at >90% / <10% median methylation. The
  feature set is therefore already coverage-conditioned — it is not a random
  sample of the methylome, which matters when interpreting variance across sites.

`RheMacAge` also carries `baboonTrimAlignExtract.sh` for the cross-species
(two-primate) clock, and the trained coefficients in `clock/ClockBetas.txt`.

---

## 📈 Unsupervised Aging Score
In the absence of raw chronological age metadata for all samples (particularly liver), we calculate an unsupervised estimated age based on global methylation dynamics:

1. Extract three global statistics per sample across all selected CpG sites: **Mean**, **Median**, and **Variability** (standard deviation).
2. Standardize these metrics ($z$-scores).
3. Compute the unsupervised **Aging Score**:
   $$\text{Aging Score} = \text{Variability}_{\text{std}} - \text{Mean}_{\text{std}} - \text{Median}_{\text{std}}$$
4. Min-max scale the resulting scores to map to a standard macaque lifespan proxy of **3 to 40 years**.

---

## 🤖 Models Compared

We generate and compare embeddings from five distinct approaches:

1. **Base Ridge**: Raw methylation percentages (baseline) with no embedding projection.
2. **Nucleotide Transformer**: 1024-dimensional sequence embeddings from HuggingFace's `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species` DNA language model, projected via matrix-vector multiplication with sample methylation levels ($M \times E$).
3. **MuLan-Methyl**: 768-dimensional sequence embeddings from the 3-mer tokenized `wenhuan/MuLan-Methyl-BERT` model, projected via $M \times E$.
4. **CpGPT Foundation Model**: A transformer/mamba architecture designed to produce sample-level representations directly from methylation beta values, genomic positions, chromosomes, and DNA embeddings.
5. **MethylBERT**: Continuous methylation embedding model utilizing a Performer-based language representation backbone to encode read-level/locus-level methylation context.


---

## 🚀 Getting Started

```bash
pip install torch transformers scikit-learn pandas numpy matplotlib seaborn scipy
jupyter notebook train_models_multimodel.ipynb
```

### Foundation model dependencies

CpGPT, MethylBERT and MuLan-Methyl are external repositories, not part of this
source tree. Clone them alongside this one and point the notebook at them:

```bash
git clone https://github.com/lucascamillomd/CpGPT.git
git clone https://github.com/methylbert/methylbert.git
git clone https://github.com/husonlab/mulan-methyl.git
```

### Data

The Cayo Santiago and ONPRC matrices are **not** in this repository — every
`*.csv` is gitignored. `prepare_data.py` rebuilds the liver subset and fetches
±100 bp `Mmul_10` context windows from the Ensembl REST API.
