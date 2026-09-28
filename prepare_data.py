import pandas as pd
import numpy as np
import json
import requests
import time
import subprocess
import os

print("Processing Liver data...")
rdata_path = 'myData.RData'
# Use R to extract the top 200 variable sites
r_code = f"""
load('{rdata_path}')
sd_vals <- apply(perc5X, 1, sd, na.rm=TRUE)
top_sites <- names(sort(sd_vals, decreasing=TRUE)[1:200])
sub_mat <- perc5X[top_sites, ]
write.csv(sub_mat, 'liver_subset.csv')
"""
with open('extract_subset.R', 'w') as f:
    f.write(r_code)

subprocess.run(['Rscript', 'extract_subset.R'])
liver_df = pd.read_csv('liver_subset.csv', index_col=0)
print(f"Extracted liver matrix: {liver_df.shape}")

print("Fetching sequences from Ensembl for macaque...")
liver_seqs = {}
for site in liver_df.index:
    parts = site.split('_')
    chrom = parts[0]
    pos = int(parts[1])
    
    # Ensembl REST API for macaque (macaca_mulatta)
    url = f"https://rest.ensembl.org/sequence/region/macaca_mulatta/{chrom}:{pos-100}..{pos+100}?coord_system_version=Mmul_10"
    try:
        r = requests.get(url, headers={"Content-Type": "application/json"})
        if r.status_code == 200:
            data = r.json()
            liver_seqs[site] = data['seq']
        else:
            print(f"Failed {site}: {r.status_code} {r.text}")
    except Exception as e:
        print(f"Error {site}: {e}")
    time.sleep(0.05)

with open('liver_selected_sequences.json', 'w') as f:
    json.dump(liver_seqs, f)

print(f"Saved {len(liver_seqs)} sequences for liver.")
