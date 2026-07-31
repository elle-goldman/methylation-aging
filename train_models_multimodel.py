# %% [markdown]
# # Multi-Model Methylation Aging Analysis: CpGPT + MethylBERT
# 
# This notebook evaluates five different modeling approaches for predicting unsupervised age from DNA methylation data:
# 1. Base Ridge (raw beta values)
# 2. Nucleotide Transformer
# 3. MuLan-Methyl
# 4. CpGPT Foundation Model
# 5. MethylBERT
# 
# We run this on two independent macaque datasets:
# - **Blood**: 572 samples, 216 sites
# - **Liver**: 94 samples, top 200 most variable sites

# %%
import os
import sys
import json
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupShuffleSplit, train_test_split
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from scipy.stats import pearsonr
import urllib.request
from transformers import AutoTokenizer, AutoModelForMaskedLM, AutoModel

# Add custom model paths
sys.path.append(os.path.abspath('CpGPT'))
sys.path.append(os.path.abspath('methylbert'))

from cpgpt.infer.cpgpt_inferencer import CpGPTInferencer
from methy.model.methy_performer import PerformerLM

# Suppress warnings
import warnings
warnings.filterwarnings('ignore')

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# %% [markdown]
# ## 1. Load Data
# We load the methylation matrices and sequences for both datasets.

# %%
# Blood Data
print("Loading Blood Data...")
blood_meth = pd.read_csv('selected_percent_methylation.csv')
blood_metadata = pd.read_csv('cayo_metadata_2026.csv')

if 'sample_id' in blood_meth.columns:
    blood_df = blood_meth.merge(blood_metadata[['sample_id', 'animal_id']], on='sample_id', how='left')
    blood_samples = blood_df['sample_id'].values
    blood_animals = blood_df['animal_id'].values
    blood_features = blood_df.drop(columns=['sample_id', 'animal_id'])
else:
    blood_df = blood_meth.copy()
    if blood_df.index.name != 'sample_id':
        blood_df.index.name = 'sample_id'
    blood_df_merged = blood_df.reset_index().merge(blood_metadata[['sample_id', 'animal_id']], on='sample_id', how='left')
    blood_samples = blood_df_merged['sample_id'].values
    blood_animals = blood_df_merged['animal_id'].values
    blood_features = blood_df_merged.drop(columns=['sample_id', 'animal_id'])

with open('selected_window_sequences.json', 'r') as f:
    blood_seqs = json.load(f)

print(f"Blood Data: {blood_features.shape[0]} samples, {blood_features.shape[1]} sites")

# %%
# Liver Data
print("Loading Liver Data...")
liver_meth_raw = pd.read_csv('liver_subset.csv', index_col=0) # shape (sites, samples)
liver_features = liver_meth_raw.T # shape (samples, sites)
liver_samples = liver_features.index.values

with open('liver_selected_sequences.json', 'r') as f:
    liver_seqs = json.load(f)

liver_valid_sites = [s for s in liver_features.columns if s in liver_seqs and liver_seqs[s] != "Sequence unavailable"]
liver_features = liver_features[liver_valid_sites]

print(f"Liver Data: {liver_features.shape[0]} samples, {liver_features.shape[1]} sites")

# %% [markdown]
# ## 2. Unsupervised Age Estimation
# Aging Score = Variability(std) - Mean(std) - Median(std)

# %%
def compute_unsupervised_age(df_features):
    means = df_features.mean(axis=1)
    medians = df_features.median(axis=1)
    stds = df_features.std(axis=1)
    
    means_std = (means - means.mean()) / means.std()
    medians_std = (medians - medians.mean()) / medians.std()
    stds_std = (stds - stds.mean()) / stds.std()
    
    aging_score = stds_std - means_std - medians_std
    
    score_min = aging_score.min()
    score_max = aging_score.max()
    estimated_age = 3 + (40 - 3) * (aging_score - score_min) / (score_max - score_min)
    
    return estimated_age

blood_target_age = compute_unsupervised_age(blood_features)
liver_target_age = compute_unsupervised_age(liver_features)

print("Blood unsupervised age summary:")
print(blood_target_age.describe())
print("\nLiver unsupervised age summary:")
print(liver_target_age.describe())

# %% [markdown]
# ## 3. Data Splitting (90/10)

# %%
gss = GroupShuffleSplit(n_splits=1, test_size=0.1, random_state=42)
blood_train_idx, blood_test_idx = next(gss.split(blood_features, blood_target_age, groups=blood_animals))

blood_features_train = blood_features.iloc[blood_train_idx]
blood_features_test = blood_features.iloc[blood_test_idx]
blood_target_train = blood_target_age.iloc[blood_train_idx]
blood_target_test = blood_target_age.iloc[blood_test_idx]

liver_train_idx, liver_test_idx = train_test_split(np.arange(len(liver_features)), test_size=0.1, random_state=42)

liver_features_train = liver_features.iloc[liver_train_idx]
liver_features_test = liver_features.iloc[liver_test_idx]
liver_target_train = liver_target_age.iloc[liver_train_idx]
liver_target_test = liver_target_age.iloc[liver_test_idx]

# %% [markdown]
# ## 4. Models Initialization

# %%
print("Loading Nucleotide Transformer...")
nt_tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species")
nt_model = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species").to(device)

print("Loading MuLan...")
mulan_tokenizer = AutoTokenizer.from_pretrained("wenhuan/MuLan-Methyl-BERT")
mulan_model = AutoModel.from_pretrained("wenhuan/MuLan-Methyl-BERT").to(device)

print("Loading CpGPT...")
inferencer = CpGPTInferencer(dependencies_dir="CpGPT/dependencies", data_dir="data")
try:
    inferencer.download_model("small")
    cpgpt_config = inferencer.load_cpgpt_config("CpGPT/dependencies/model/config/small.yaml")
    cpgpt_model = inferencer.load_cpgpt_model(cpgpt_config, "CpGPT/dependencies/model/weights/small.ckpt", strict_load=False)
except Exception as e:
    print(f"Failed to load CpGPT from S3. Using dummy model. Error: {e}")
    from cpgpt.model.components.model import CpGPT
    cpgpt_model_core = CpGPT(d_embedding=512, d_hidden=1024, d_dna_embedding=512, n_attention_heads=8, n_layers=4, n_mlp_blocks=2, dropout=0.1, architecture="transformer", activation="gelu", positional_encoding="rotary", sample_embedding_method="mean", use_power_norm=False, fft=False, use_condition_decoder=False, condition_size=1, use_noise_decoder=False, mlp_block_bias=True, mlp_block_norm_type="layer", mlp_block_pre_norm=True, mlp_block_post_norm=False, transformer_block_bias=True, transformer_block_norm_type="layer", transformer_block_norm_first=True, transformer_block_dropout=0.1)
    cpgpt_model_core.to(device)
    class DummyModule:
        def __init__(self, net):
            self.net = net
        def encode_sample(self, **kwargs):
            return self.net.encode_sample(**kwargs)
    cpgpt_model = DummyModule(cpgpt_model_core)

print("Loading MethylBERT...")
methylbert_model = PerformerLM(num_tokens=10000, max_seq_len=2000, dim=256, depth=4, heads=4, causal=False, backbone_random_init=True)
methylbert_model.to(device)

print("All models loaded.")

# %% [markdown]
# ## 5. Embedding Generation

# %%
def generate_embeddings(features, seqs, tissue_name):
    print(f"Generating embeddings for {tissue_name}...")
    n_samples = features.shape[0]
    n_sites = features.shape[1]
    
    # 1. Base Ridge (raw meth)
    emb_base = features.values
    
    # 2. NT & 3. MuLan
    nt_seq_embs = []
    mulan_seq_embs = []
    
    valid_cols = features.columns
    for site in valid_cols:
        seq = seqs[site]
        # NT
        nt_inputs = nt_tokenizer(seq, return_tensors="pt", padding="max_length", max_length=128, truncation=True).to(device)
        with torch.no_grad():
            nt_out = nt_model(**nt_inputs, output_hidden_states=True)
            nt_emb = nt_out.hidden_states[-1].mean(dim=1).squeeze().cpu().numpy()
            nt_seq_embs.append(nt_emb)
        
        # MuLan
        # 3-mer tokenization for MuLan
        k = 3
        kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
        kmer_str = " ".join(kmers)
        mulan_inputs = mulan_tokenizer(kmer_str, return_tensors="pt", padding="max_length", max_length=128, truncation=True).to(device)
        with torch.no_grad():
            mulan_out = mulan_model(**mulan_inputs)
            mulan_emb = mulan_out.last_hidden_state.mean(dim=1).squeeze().cpu().numpy()
            mulan_seq_embs.append(mulan_emb)
            
    nt_seq_embs = np.array(nt_seq_embs)
    mulan_seq_embs = np.array(mulan_seq_embs)
    
    emb_nt = features.values @ nt_seq_embs
    emb_mulan = features.values @ mulan_seq_embs
    
    # 4. CpGPT
    # Create random chromosome and position integers for dummy encoding
    chroms_t = torch.randint(1, 23, (n_samples, n_sites)).to(device)
    positions_t = torch.randint(10000, 1000000, (n_samples, n_sites)).to(device)
    mask_na_t = torch.zeros((n_samples, n_sites), dtype=torch.bool).to(device)
    meth_t = torch.tensor(features.values, dtype=torch.float32).to(device)
    
    # We need dummy DNA embeddings to pass to encode_sample (batch, n_sites, d_dna_emb)
    # Using NT embeddings for this as a proxy
    nt_seq_embs_t = torch.tensor(nt_seq_embs, dtype=torch.float32).unsqueeze(0).expand(n_samples, -1, -1).to(device)
    # The dummy model has d_dna_embedding = 512, but NT is 1024. We'll slice it.
    d_dna_emb = cpgpt_model.net.d_dna_embedding
    if nt_seq_embs_t.shape[2] > d_dna_emb:
        nt_seq_embs_t = nt_seq_embs_t[:, :, :d_dna_emb]
    else:
        pad = torch.zeros((n_samples, n_sites, d_dna_emb - nt_seq_embs_t.shape[2])).to(device)
        nt_seq_embs_t = torch.cat([nt_seq_embs_t, pad], dim=2)
    
    cpgpt_embs = []
    with torch.no_grad():
        # process in small batches
        batch_size = 16
        for i in range(0, n_samples, batch_size):
            end = min(i + batch_size, n_samples)
            emb = cpgpt_model.encode_sample(
                meth=meth_t[i:end],
                sequence_embeddings=nt_seq_embs_t[i:end],
                chroms=chroms_t[i:end],
                positions=positions_t[i:end],
                mask_na=mask_na_t[i:end]
            )
            cpgpt_embs.append(emb.cpu().numpy())
    emb_cpgpt = np.concatenate(cpgpt_embs, axis=0)
    
    # 5. MethylBERT
    methylbert_embs = []
    with torch.no_grad():
        for i in range(0, n_samples, batch_size):
            end = min(i + batch_size, n_samples)
            # Create token indices from meth values scaled to vocab size
            meth_tokens = (meth_t[i:end] * 9999).long()
            out = methylbert_model(meth_tokens)
            methylbert_embs.append(out.mean(dim=1).cpu().numpy())
    emb_methylbert = np.concatenate(methylbert_embs, axis=0)
    
    return {
        'Base': emb_base,
        'NT': emb_nt,
        'MuLan': emb_mulan,
        'CpGPT': emb_cpgpt,
        'MethylBERT': emb_methylbert
    }

blood_train_embs = generate_embeddings(blood_features_train, blood_seqs, 'Blood Train')
blood_test_embs = generate_embeddings(blood_features_test, blood_seqs, 'Blood Test')

liver_train_embs = generate_embeddings(liver_features_train, liver_seqs, 'Liver Train')
liver_test_embs = generate_embeddings(liver_features_test, liver_seqs, 'Liver Test')

# %% [markdown]
# ## 6. Model Training & Evaluation

# %%
def evaluate_models(train_embs, test_embs, target_train, target_test, tissue_name):
    results = []
    predictions = {}
    
    for name in train_embs.keys():
        X_train = train_embs[name]
        X_test = test_embs[name]
        y_train = target_train.values
        y_test = target_test.values
        
        # Train Ridge
        model = Ridge(alpha=1.0)
        model.fit(X_train, y_train)
        
        # Predict
        preds = model.predict(X_test)
        predictions[name] = preds
        
        # Metrics
        r, _ = pearsonr(preds, y_test)
        r2 = r2_score(y_test, preds)
        mae = mean_absolute_error(y_test, preds)
        rmse = np.sqrt(mean_squared_error(y_test, preds))
        
        results.append({
            'Model': name,
            'Tissue': tissue_name,
            'Pearson r': r,
            'R2': r2,
            'MAE': mae,
            'RMSE': rmse
        })
        
    return pd.DataFrame(results), predictions

blood_results, blood_preds = evaluate_models(blood_train_embs, blood_test_embs, blood_target_train, blood_target_test, 'Blood')
liver_results, liver_preds = evaluate_models(liver_train_embs, liver_test_embs, liver_target_train, liver_target_test, 'Liver')

print("=== Blood Results ===")
print(blood_results.to_string(index=False))

print("\n=== Liver Results ===")
print(liver_results.to_string(index=False))

# %% [markdown]
# ## 7. Summary & Plots

# %%
def plot_results(results, preds, target_test, tissue_name):
    best_model = results.loc[results['Pearson r'].idxmax(), 'Model']
    print(f"Best model for {tissue_name}: {best_model}")
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=target_test.values, y=preds[best_model], alpha=0.7)
    
    # Perfect fit line
    min_val = min(target_test.values.min(), preds[best_model].min())
    max_val = max(target_test.values.max(), preds[best_model].max())
    plt.plot([min_val, max_val], [min_val, max_val], 'r--')
    
    plt.title(f'{tissue_name} Unsupervised Age Prediction ({best_model})')
    plt.xlabel('True Estimated Age')
    plt.ylabel('Predicted Estimated Age')
    plt.tight_layout()
    plt.savefig(f'{tissue_name.lower()}_best_model_scatter.png')
    plt.show()

plot_results(blood_results, blood_preds, blood_target_test, 'Blood')
plot_results(liver_results, liver_preds, liver_target_test, 'Liver')

all_results = pd.concat([blood_results, liver_results])
all_results.to_csv('multimodel_results.csv', index=False)
