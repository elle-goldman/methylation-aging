# 🧬 Methylation Aging: Multi-Model Epigenetic Clock Analysis

This repository implements a comparative framework evaluating multiple DNA methylation foundation models and sequence-based models for predicting unsupervised biological age across independent rhesus macaque (*Macaca mulatta*) tissues.

---

## 📋 Table of Contents
1. [Project Overview](#-project-overview)
2. [Datasets](#-datasets)
3. [Unsupervised Aging Score](#-unsupervised-aging-score)
4. [Models Compared](#-models-compared)
5. [Project Structure](#-project-structure)
6. [Dependencies & Credits](#-dependencies--credits)
7. [Getting Started](#-getting-started)

---

## 🧬 Project Overview
Epigenetic clocks traditionally map DNA methylation (DNAm) levels at specific CpG sites to chronological age. This project compares standard regression against representations derived from DNA language models (LLMs) and methylation transformer models. 

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

## 📂 Project Structure
```
├── CpGPT/                      # Submodule: CpGPT model repository
├── methylbert/                 # Submodule: MethylBERT model repository
├── mulan-methyl/               # Submodule: MuLan-Methyl model repository
├── selected_percent_methylation.csv  # Blood methylation levels (572 x 216)
├── selected_window_sequences.json    # DNA sequence windows for blood CpG sites
├── liver_subset.csv            # Pre-extracted top 200 variable liver CpG sites
├── liver_selected_sequences.json     # DNA sequence windows for liver CpG sites
├── cayo_metadata_2026.csv      # Metadata for Cayo blood dataset
├── prepare_data.py             # Script to extract liver subset & fetch sequences
├── train_models_multimodel.ipynb # Primary analysis & evaluation notebook
├── train_models_multimodel.py    # Python conversion of the analysis script
└── README.md                   # This file
```

---

## 🤝 Models Compared

*   **CpGPT**:
    *   *Paper*: "CpGPT: A Foundation Model for DNA Methylation"
    *   *Repository*: [lucascamillomd/CpGPT](https://github.com/lucascamillomd/CpGPT)
*   **MethylBERT**:
    *   *Paper*: "methylBERT: A transformer-based model for read-level DNA methylation patterns"
    *   *Repository*: [methylbert/methylbert](https://github.com/methylbert/methylbert)
*   **MuLan-Methyl**:
    *   *Paper*: "MuLan-Methyl: Multi-Task Language Representation Fine-Tuning for DNA Methylation Prediction"
    *   *Repository*: [husonlab/mulan-methyl](https://github.com/husonlab/mulan-methyl)

---

## 🚀 Getting Started

### 1. Clone the Repository
Because this repository tracks submodules, you must clone recursively to fetch the dependencies:
```bash
git clone --recursive <repository-url>
cd methylation-aging
```

If you already cloned the repository without submodules, run:
```bash
git submodule update --init --recursive
```

### 2. Install Environment & Requirements
Ensure you have the required Python environments. The notebook relies on PyTorch, HuggingFace Transformers, scikit-learn, and other common libraries.
```bash
pip install torch transformers scikit-learn pandas numpy matplotlib seaborn scipy
```

### 3. Run the Analysis
Open the Jupyter notebook and run the cells sequentially to extract embeddings, run Ridge regression cross-validation, and compare performance tables and scatter plots:
```bash
jupyter notebook train_models_multimodel.ipynb
```
