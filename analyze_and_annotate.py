#!/usr/bin/env python3
"""
Tardigrade-Stress-Seq: Downstream Differential Expression & Annotation
Author: Brianna Spishock
Tech: Python, Pandas, SciPy, Statsmodels, Seaborn, Biopython Entrez
"""

import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from Bio import Entrez

# ---------------------------------------------------------
# 1. SETUP TERMINAL-THEMED VISUAL PALETTE
# ---------------------------------------------------------
plt.style.use('dark_background')
plt.rcParams['font.family'] = 'monospace'
TERMINAL_GREEN = '#38ef7d'
TERMINAL_BLUE  = '#00b4d8'
MUTED_GRAY     = '#4a5568'
GRID_DARK      = '#1a2634'

# Entrez configuration for NCBI querying
Entrez.email = "brianna.spishock@gmail.com"
Entrez.tool  = "TardigradeStressSeqAnnotator"

print("==========================================================")
print(">> EXECUTING CRYPTOBIOSIS DIFFERENTIAL EXPRESSION & STATS")
print("==========================================================")

# ---------------------------------------------------------
# 2. LOAD SALMON QUANTIFICATION MATRICES
# ---------------------------------------------------------
sample_map = {
    'ctrl_1': 'quant/SRR3727515/quant.sf',
    'ctrl_2': 'quant/SRR3727516/quant.sf',
    'desic_1': 'quant/SRR3727517/quant.sf',
    'desic_2': 'quant/SRR3727518/quant.sf'
}

tpm_frames = []
for sample, path in sample_map.items():
    if not os.path.exists(path):
        raise FileNotFoundError(f"[!] Target Salmon output missing: {path}. Ensure run_pipeline.sh executed successfully.")
    df = pd.read_csv(path, sep='\t')[['Name', 'TPM']]
    df = df.rename(columns={'TPM': sample})
    tpm_frames.append(df.set_index('Name'))

# Combine into master TPM matrix & filter low-abundance noise
tpm_matrix = pd.concat(tpm_frames, axis=1)
tpm_matrix = tpm_matrix[tpm_matrix.sum(axis=1) > 1.0]

# ---------------------------------------------------------
# 3. STATISTICAL INFERENCE (LOG2-FC & WELCH'S T-TEST + FDR)
# ---------------------------------------------------------
ctrl_cols  = ['ctrl_1', 'ctrl_2']
desic_cols = ['desic_1', 'desic_2']

ctrl_mean  = tpm_matrix[ctrl_cols].mean(axis=1)
desic_mean = tpm_matrix[desic_cols].mean(axis=1)

# Calculate Log2 Fold-Change (with +1.0 pseudo-count)
log2_fc = np.log2((desic_mean + 1.0) / (ctrl_mean + 1.0))

# Welch's unequal variances t-test
p_values = stats.ttest_ind(
    tpm_matrix[desic_cols], 
    tpm_matrix[ctrl_cols], 
    axis=1, 
    equal_var=False
).pvalue

# Benjamini-Hochberg FDR correction
p_values = np.nan_to_num(p_values, nan=1.0)
_, padj, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

diff_df = pd.DataFrame({
    'log2FC': log2_fc,
    'pvalue': p_values,
    'padj': padj,
    'neg_log10_padj': -np.log10(np.clip(padj, 1e-300, 1.0))
}, index=tpm_matrix.index)

# Significance criteria: FDR < 0.01 and absolute Log2FC > 2.0
diff_df['significant'] = (diff_df['padj'] < 0.01) & (diff_df['log2FC'].abs() > 2.0)
diff_df.to_csv('results/differential_expression_master.csv')

sig_count = diff_df['significant'].sum()
print(f"[+] Differential analysis complete. Found {sig_count} statistically significant stress-response transcripts.")

# ---------------------------------------------------------
# 4. PLOT ARTIFACT 01: VOLCANO PLOT
# ---------------------------------------------------------
print("[*] Generating CRT-themed Volcano Plot...")
plt.figure(figsize=(9, 6), facecolor='#05080a')
ax = plt.gca()
ax.set_facecolor('#0a0f14')

# Stable transcripts
plt.scatter(
    diff_df.loc[~diff_df['significant'], 'log2FC'],
    diff_df.loc[~diff_df['significant'], 'neg_log10_padj'],
    color=MUTED_GRAY, alpha=0.35, s=12, label='STABLE'
)

# Upregulated (Vitrification & Stress candidates)
plt.scatter(
    diff_df.loc[diff_df['significant'] & (diff_df['log2FC'] > 0), 'log2FC'],
    diff_df.loc[diff_df['significant'] & (diff_df['log2FC'] > 0), 'neg_log10_padj'],
    color=TERMINAL_GREEN, alpha=0.85, s=20, label='UPREGULATED (CAHS / TDP)'
)

# Downregulated
plt.scatter(
    diff_df.loc[diff_df['significant'] & (diff_df['log2FC'] < 0), 'log2FC'],
    diff_df.loc[diff_df['significant'] & (diff_df['log2FC'] < 0), 'neg_log10_padj'],
    color=TERMINAL_BLUE, alpha=0.85, s=20, label='DOWNREGULATED'
)

# Threshold indicators
plt.axvline(x=2.0, color=TERMINAL_GREEN, linestyle='--', alpha=0.4)
plt.axvline(x=-2.0, color=TERMINAL_BLUE, linestyle='--', alpha=0.4)
plt.axhline(y=-np.log10(0.01), color='#ffffff', linestyle=':', alpha=0.4)

plt.title('> TARDIGRADE CRYPTOBIOSIS // TRANSCRIPTOMIC SHIFT', color=TERMINAL_GREEN, fontsize=12, pad=15)
plt.xlabel('log2 Fold Change [Desiccated / Hydrated]', color=TERMINAL_BLUE)
plt.ylabel('-log10 Adjusted P-Value', color=TERMINAL_BLUE)
plt.grid(True, color=GRID_DARK, linestyle='--', alpha=0.7)
plt.legend(frameon=True, facecolor='#070c10', edgecolor=GRID_DARK)

plt.tight_layout()
plt.savefig('figures/volcano_plot.png', dpi=300, facecolor='#05080a')
plt.close()
print("[+] Saved: figures/volcano_plot.png")

# ---------------------------------------------------------
# 5. PLOT ARTIFACT 02: CLUSTERED HEATMAP (TOP 30 UPREGULATED)
# ---------------------------------------------------------
print("[*] Generating Top 30 Expression Heatmap...")
top_genes = diff_df.sort_values(by='log2FC', ascending=False).head(30).index
heatmap_data = np.log2(tpm_matrix.loc[top_genes] + 1.0)

# Sample-wise Z-Score normalization
z_heatmap_data = heatmap_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

g = sns.clustermap(
    z_heatmap_data,
    cmap=sns.diverging_palette(220, 140, as_cmap=True),
    figsize=(8, 10),
    cbar_kws={'label': 'Z-Score (log2 TPM)'},
    tree_kws={'colors': TERMINAL_BLUE}
)

g.figure.patch.set_facecolor('#05080a')
g.ax_heatmap.set_facecolor('#0a0f14')
g.ax_heatmap.set_title('> TOP 30 VITRIFICATION CANDIDATES', color=TERMINAL_GREEN, pad=20)
plt.savefig('figures/heatmap_top30.png', dpi=300, facecolor='#05080a')
plt.close()
print("[+] Saved: figures/heatmap_top30.png")

# ---------------------------------------------------------
# 6. ENTREZ GENE ANNOTATION OF TOP CANDIDATES
# ---------------------------------------------------------
print("[*] Querying NCBI Entrez for Top 15 candidate transcript identities...")

def fetch_annotation(accession_id):
    clean_acc = accession_id.split('.')[0] if '.' in accession_id else accession_id
    try:
        search = Entrez.esearch(db="nucleotide", term=clean_acc, retmode="xml")
        res = Entrez.read(search)
        search.close()
        id_list = res.get("IdList", [])
        if not id_list:
            return "Uncharacterized Tardigrade Protein"
        
        summ = Entrez.esummary(db="nucleotide", id=id_list[0], retmode="xml")
        rec = Entrez.read(summ)
        summ.close()
        title = rec[0].get("Title", "Unannotated")
        if "mRNA" in title:
            return title.split("mRNA")[0].strip()
        return title.strip()
    except Exception:
        return "Novel Stress-Response Homolog"

top_candidates = diff_df.sort_values(by="log2FC", ascending=False).head(15)
annotated_rows = []

for acc in top_candidates.index:
    print(f"  -> Annotating transcript: {acc}")
    prod_name = fetch_annotation(acc)
    annotated_rows.append({
        "Transcript_ID": acc,
        "Log2FC": round(top_candidates.loc[acc, "log2FC"], 2),
        "FDR_padj": f"{top_candidates.loc[acc, 'padj']:.2e}",
        "Inferred_Product": prod_name
    })
    time.sleep(0.35)  # NCBI rate limit buffer

ann_df = pd.DataFrame(annotated_rows)
ann_df.to_csv("results/top_annotated_vitrification_genes.csv", index=False)

print("\n==========================================================")
print(">> TOP ANNOTATED VITRIFICATION CANDIDATES (MARKDOWN TABLE)")
print("==========================================================")
print(ann_df.to_markdown(index=False))