import os
import json
import pandas as pd
import requests
from pyliftover import LiftOver

print("Initializing LiftOver objects...")
scratch_dir = '/Users/hellyr/.gemini/antigravity-ide/brain/5edd0dc4-3f77-4671-b946-3300bc9e5423/scratch/'
dest_10 = os.path.join(scratch_dir, 'rheMac10ToHg38.over.chain.gz')
dest_8 = os.path.join(scratch_dir, 'rheMac8ToHg38.over.chain.gz')

lo_10 = LiftOver(dest_10)
lo_8 = LiftOver(dest_8)
print("LiftOver loaded successfully!")

# Load top correlations
top_cors = pd.read_csv('macaque_top_correlations.csv')

results = []

for idx, row in top_cors.iterrows():
    cpg = row['cpg']
    correlation = row['correlation']
    
    parts = cpg.split('_')
    chrom = 'chr' + parts[0]
    pos = int(parts[1])
    
    # Try rheMac10 first, then rheMac8
    mapped_coord = lo_10.convert_coordinate(chrom, pos)
    assembly_used = 'rheMac10'
    if not mapped_coord:
        mapped_coord = lo_8.convert_coordinate(chrom, pos)
        assembly_used = 'rheMac8'
        
    if mapped_coord:
        human_chrom = mapped_coord[0][0]
        human_pos = mapped_coord[0][1]
        human_strand = mapped_coord[0][2]
        
        # Ensembl REST API query to find overlapping genes
        # We query a region of +/- 5000 bp around the CpG site to find nearby genes
        region_str = f"{human_chrom.replace('chr', '')}:{human_pos-5000}-{human_pos+5000}"
        url = f"https://rest.ensembl.org/overlap/region/human/{region_str}?feature=gene"
        
        genes = []
        try:
            r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=10)
            if r.status_code == 200:
                gene_data = r.json()
                for g in gene_data:
                    gene_name = g.get('external_name', g.get('gene_id', 'Unknown'))
                    gene_id = g.get('id', '')
                    biotype = g.get('biotype', '')
                    
                    # Query GO terms for this gene ID
                    go_url = f"https://rest.ensembl.org/xrefs/id/{gene_id}?all_levels=1"
                    go_terms = []
                    r_go = requests.get(go_url, headers={"Content-Type": "application/json"}, timeout=10)
                    if r_go.status_code == 200:
                        go_data = r_go.json()
                        for term in go_data:
                            if term.get('dbname') == 'GO':
                                go_terms.append({
                                    'go_id': term.get('primary_id', ''),
                                    'go_name': term.get('description', '')
                                })
                    
                    genes.append({
                        'gene_name': gene_name,
                        'gene_id': gene_id,
                        'biotype': biotype,
                        'go_terms': go_terms
                    })
        except Exception as e:
            print(f"Error querying Ensembl for {region_str}: {e}")
            
        results.append({
            'macaque_cpg': cpg,
            'correlation': correlation,
            'assembly': assembly_used,
            'human_coords': f"{human_chrom}:{human_pos}",
            'genes': genes
        })
    else:
        results.append({
            'macaque_cpg': cpg,
            'correlation': correlation,
            'assembly': 'Failed',
            'human_coords': 'Failed to liftOver',
            'genes': []
        })

# Save JSON results
with open('cpg_liftover_results.json', 'w') as f:
    json.dump(results, f, indent=2)

# Generate Markdown Summary
md_content = """# LiftOver and GO Annotation of Top Age-Predictive CpG Sites

This report details the conversion of the top 20 rhesus macaque CpG sites (predictive of liver aging and health) to the human genome reference (GRCh38/hg38) and their biological annotation.

## Summary Table:
#$ Contents
| Macaque CpG | Correlation | Assembly | Human Coords (hg38) | Nearby/Overlapping Genes | Top GO Terms |
|---|---|---|---|---|---|
"""

for r in results:
    cpg = r['macaque_cpg']
    corr = f"{r['correlation']:.3f}"
    asm = r['assembly']
    h_coords = r['human_coords']
    
    if r['genes']:
        gene_names = ", ".join([g['gene_name'] for g in r['genes']])
        
        # Collect GO terms
        all_go = []
        for g in r['genes']:
            for term in g['go_terms']:
                if term['go_name'] and term['go_name'] not in all_go:
                    all_go.append(term['go_name'])
        top_go = "; ".join(all_go[:3]) if all_go else "None annotated"
    else:
        gene_names = "None"
        top_go = "None"
        
    md_content += f"| {cpg} | {corr} | {asm} | {h_coords} | {gene_names} | {top_go} |\n"

with open('/Users/hellyr/Desktop/Projects/Methylation Aging/cpg_liftover_summary.md', 'w') as f:
    f.write(md_content)

print("Mapping and annotation finished successfully!")
