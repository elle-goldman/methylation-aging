"""
Rerun Cayo Santiago macaque aging analysis with CpG density incorporated
as an additional feature alongside CpGPT (Nucleotide Transformer) embeddings.

CpG density = count("CG") / len(sequence) for each window's 2001bp context.
Per-sample CpG density score = M @ cpg_density_vector (weighted sum of window densities).
This is appended as extra feature(s) to the 1024-dim CpGPT embeddings.
"""

import os
import json
import pandas as pd
import numpy as np
import torch

# Inject missing find_pruneable_heads_and_indices into transformers.pytorch_utils for v5 compatibility
import transformers.pytorch_utils
def find_pruneable_heads_and_indices(heads, n_heads, head_size, already_pruned_heads):
    mask = torch.ones(n_heads, head_size)
    heads = set(heads) - already_pruned_heads
    for head in heads:
        mask[head] = 0
    mask = mask.view(-1)
    index = mask.nonzero().squeeze(-1)
    return heads, index
transformers.pytorch_utils.find_pruneable_heads_and_indices = find_pruneable_heads_and_indices

from transformers import AutoTokenizer, AutoModelForMaskedLM, AutoConfig, PreTrainedModel, AutoModel
PreTrainedModel.all_tied_weights_keys = {}

def get_head_mask(self, head_mask, num_hidden_layers, is_attention_chunked=False):
    if head_mask is None:
        return [None] * num_hidden_layers
    return head_mask
PreTrainedModel.get_head_mask = get_head_mask

from sklearn.model_selection import GroupKFold
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from scipy.stats import pearsonr

# ============================================================
# File paths
# ============================================================
sequences_path = '/Users/hellyr/Desktop/Projects/Methylation Aging/selected_window_sequences.json'
meth_path = '/Users/hellyr/Desktop/Projects/Methylation Aging/selected_percent_methylation.csv'
metadata_path = '/Users/hellyr/Desktop/Projects/Methylation Aging/cayo_metadata_2026.csv'
output_results_csv = '/Users/hellyr/Desktop/Projects/Methylation Aging/cayo_model_predictions_v2.csv'

# ============================================================
# 1. Load sequences and compute CpG density per window
# ============================================================
print("Loading sequences...")
with open(sequences_path, 'r') as f:
    window_data = json.load(f)

selected_windows_in_order = [row['macaque_window'] for row in window_data]
cpgpt_seqs = [row['cpgpt_seq'] for row in window_data]
mulan_seqs = [row['mulan_seq'] for row in window_data]

# Compute CpG density for each window's 2001bp CpGPT context
cpg_densities = []
for seq in cpgpt_seqs:
    if not seq or len(seq) < 50:
        cpg_densities.append(0.0)
    else:
        cpg_count = seq.upper().count('CG')
        cpg_densities.append(cpg_count / len(seq))
cpg_densities = np.array(cpg_densities)

print(f"Loaded {len(window_data)} window sequences.")
valid_mask = cpg_densities > 0
print(f"  Valid sequences: {valid_mask.sum()}, missing/short: {(~valid_mask).sum()}")
print(f"  CpG density range: {cpg_densities[valid_mask].min():.4f} - {cpg_densities[valid_mask].max():.4f}")
print(f"  CpG density mean:  {cpg_densities[valid_mask].mean():.4f} ± {cpg_densities[valid_mask].std():.4f}")

# ============================================================
# 2. Compute CpGPT (Nucleotide Transformer) Embeddings
# ============================================================
device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))
print(f"\nUsing device: {device}")

print("\nLoading CpGPT (Nucleotide Transformer)...")
cpgpt_model_name = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
cpgpt_tokenizer = AutoTokenizer.from_pretrained(cpgpt_model_name, trust_remote_code=True)
cpgpt_config = AutoConfig.from_pretrained(cpgpt_model_name, trust_remote_code=True)

def patch_getattr(self, name):
    defaults = {
        'rope_theta': 10000.0,
        'is_decoder': False,
        'add_cross_attention': False,
        'is_encoder_decoder': False,
        'chunk_size_feed_forward': 0,
        'output_attentions': False,
        'output_hidden_states': False,
        'use_return_dict': True,
        'pruned_heads': {},
        'tie_word_embeddings': True
    }
    if name in defaults:
        return defaults[name]
    raise AttributeError(f"'EsmConfig' object has no attribute '{name}'")

cpgpt_config.__class__.__getattr__ = patch_getattr

cpgpt_model = AutoModelForMaskedLM.from_pretrained(cpgpt_model_name, config=cpgpt_config, trust_remote_code=True).to(device)
cpgpt_model.eval()

print("Generating CpGPT embeddings...")
cpgpt_embs = []
for idx, seq in enumerate(cpgpt_seqs):
    if not seq or len(seq) < 50:
        cpgpt_embs.append(np.zeros(1024))
        continue
    inputs = cpgpt_tokenizer(seq, return_tensors="pt", truncation=True, max_length=2048)
    inputs = {k: v.to(device) for k, v in inputs.items()}
    with torch.no_grad():
        outputs = cpgpt_model(**inputs, output_hidden_states=True)
        hidden = outputs.hidden_states[-1]
        mean_emb = hidden.mean(dim=1).squeeze().cpu().numpy()
        cpgpt_embs.append(mean_emb)

E_cpgpt = np.array(cpgpt_embs)  # Shape: (216, 1024)
print(f"CpGPT sequence embedding matrix shape: {E_cpgpt.shape}")

# ============================================================
# 3. Compute MuLan-Methyl Embeddings
# ============================================================
print("\nLoading MuLan-Methyl BERT...")
mulan_model_name = "wenhuan/MuLan-Methyl-BERT"
mulan_tokenizer = AutoTokenizer.from_pretrained(mulan_model_name, trust_remote_code=True)
mulan_model = AutoModel.from_pretrained(mulan_model_name, trust_remote_code=True).to(device)
mulan_model.eval()

print("Generating MuLan-Methyl embeddings...")
mulan_embs = []
for idx, seq in enumerate(mulan_seqs):
    if not seq or len(seq) < 10:
        mulan_embs.append(np.zeros(768))
        continue
    inputs = mulan_tokenizer(seq, return_tensors="pt", truncation=True, max_length=128)
    inputs = {k: v.to(device) for k, v in inputs.items()}
    with torch.no_grad():
        outputs = mulan_model(**inputs, output_hidden_states=True)
        hidden = outputs.hidden_states[-1]
        mean_emb = hidden.mean(dim=1).squeeze().cpu().numpy()
        mulan_embs.append(mean_emb)

E_mulan = np.array(mulan_embs)  # Shape: (216, 768)
print(f"MuLan sequence embedding matrix shape: {E_mulan.shape}")

# ============================================================
# 4. Load Methylation and Match Samples
# ============================================================
print("\nLoading methylation matrix...")
meth_df = pd.read_csv(meth_path, index_col=0)
meth_df = meth_df[selected_windows_in_order]
M = meth_df.values  # Shape: (572, 216)
print(f"Methylation matrix shape: {M.shape}")

# Original sample embeddings (baseline, no CpG density)
S_cpgpt_orig = np.dot(M, E_cpgpt)     # (572, 1024)
S_mulan = np.dot(M, E_mulan)           # (572, 768)

# CpG-density-as-feature: compute per-sample CpG density score
# Each sample gets a weighted sum of CpG densities, weighted by its methylation values
cpg_density_score = np.dot(M, cpg_densities).reshape(-1, 1)  # (572, 1)
print(f"Per-sample CpG density score shape: {cpg_density_score.shape}")
print(f"  Range: {cpg_density_score.min():.4f} - {cpg_density_score.max():.4f}")
print(f"  Mean:  {cpg_density_score.mean():.4f} ± {cpg_density_score.std():.4f}")

# Append CpG density score to CpGPT embeddings
S_cpgpt_density = np.hstack([S_cpgpt_orig, cpg_density_score])  # (572, 1025)
print(f"CpGPT + CpG density feature matrix shape: {S_cpgpt_density.shape}")

# Load metadata
meta_df = pd.read_csv(metadata_path)
meta_df.set_index('Unnamed: 0', inplace=True)
y = meta_df.loc[meth_df.index, 'age_at_sample'].values
groups = meta_df.loc[meth_df.index, 'animal_id'].values

print(f"\nMerged {len(y)} samples with true chronological age.")

# ============================================================
# 5. Grouped 5-Fold Cross-Validation
# ============================================================
def run_cv(X_features, y_target, groups_list, model_name=""):
    gkf = GroupKFold(n_splits=5)
    preds_out = np.zeros(len(y_target))

    for fold, (train_idx, test_idx) in enumerate(gkf.split(X_features, y_target, groups=groups_list)):
        X_train, X_test = X_features[train_idx], X_features[test_idx]
        y_train, y_test = y_target[train_idx], y_target[test_idx]

        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)

        model = Ridge(alpha=10.0)
        model.fit(X_train_scaled, y_train)
        preds_out[test_idx] = model.predict(X_test_scaled)

    # Calculate overall metrics
    r_corr, _ = pearsonr(y_target, preds_out)
    r2 = r2_score(y_target, preds_out)
    mae = mean_absolute_error(y_target, preds_out)
    rmse = np.sqrt(mean_squared_error(y_target, preds_out))

    print(f"  {model_name}: r={r_corr:.4f}, R²={r2:.4f}, MAE={mae:.4f}, RMSE={rmse:.4f}")
    return preds_out, r_corr, r2, mae, rmse


print("\n" + "="*70)
print("Running Grouped 5-Fold Cross-Validation")
print("="*70)

print("\n--- Base Elastic Net (methylation features only) ---")
en_cv_preds, r_en, r2_en, mae_en, rmse_en = run_cv(M, y, groups, "ElasticNet Base")

print("\n--- CpGPT Original (no CpG density) ---")
cpgpt_orig_preds, r_orig, r2_orig, mae_orig, rmse_orig = run_cv(S_cpgpt_orig, y, groups, "CpGPT Original")

print("\n--- CpGPT + CpG Density Feature ---")
cpgpt_density_preds, r_dens, r2_dens, mae_dens, rmse_dens = run_cv(S_cpgpt_density, y, groups, "CpGPT + CpG Density")

print("\n--- MuLan-Methyl ---")
mulan_cv_preds, r_mulan, r2_mulan, mae_mulan, rmse_mulan = run_cv(S_mulan, y, groups, "MuLan-Methyl")

# ============================================================
# 6. Save Results
# ============================================================
output_df = pd.DataFrame({
    'sample': meth_df.index,
    'true_age': y,
    'animal_id': groups,
    'en_predicted_age': en_cv_preds,
    'cpgpt_predicted_age': cpgpt_orig_preds,
    'cpgpt_density_predicted_age': cpgpt_density_preds,
    'mulan_predicted_age': mulan_cv_preds
})
output_df.to_csv(output_results_csv, index=False)
print(f"\nPredictions saved to {output_results_csv}")

# ============================================================
# 7. Summary Comparison Table
# ============================================================
print("\n" + "="*70)
print("MODEL PERFORMANCE COMPARISON")
print("="*70)

summary_data = {
    "Model": [
        "Base Elastic Net (216 meth features)",
        "CpGPT Original (1024-dim, no density)",
        "CpGPT + CpG Density Feature (1025-dim)",
        "MuLan-Methyl (768-dim)"
    ],
    "Pearson r": [r_en, r_orig, r_dens, r_mulan],
    "R²": [r2_en, r2_orig, r2_dens, r2_mulan],
    "MAE (yrs)": [mae_en, mae_orig, mae_dens, mae_mulan],
    "RMSE (yrs)": [rmse_en, rmse_orig, rmse_dens, rmse_mulan]
}
summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))

# Highlight improvement
delta_r = r_dens - r_orig
delta_mae = mae_dens - mae_orig
print(f"\nCpG Density Feature Impact (vs original CpGPT):")
print(f"  Δ Pearson r:  {delta_r:+.4f}")
print(f"  Δ MAE:        {delta_mae:+.4f} years")
print(f"  Δ R²:         {r2_dens - r2_orig:+.4f}")
