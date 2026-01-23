import anndata as ad
import scanpy as sc

carmdr = sc.read_10x_h5('hao_hcc/cardrm/filtered_feature_bc_matrix.h5')
carmdr.var_names_make_unique()  # Ensure unique gene names

import scanpy as sc
import os

# Optional: Confirm current directory
print(os.getcwd())  # Should show /storage/riosugimuralab/alexto/B1_objects/B1_Sx/hao_hcc/carm

carm = sc.read_10x_mtx('hao_hcc/carm', var_names='gene_symbols', cache=True)
carm.var_names_make_unique()



# Assume carm and carmdr are loaded and preprocessed similarly (QC, normalized)
carm.obs['sample'] = 'carm'
carmdr.obs['sample'] = 'carmdr'

adata = ad.concat([carm, carmdr], 
                  label='sample', 
                  keys=['carm', 'carmdr'], 
                  join='outer')  # Keeps union of genes
adata.var_names_make_unique()


import harmonypy as hm
import numpy as np

# Use harmonypy directly for more control
data_mat = adata.obsm['X_pca'].astype(np.float64).T  # Harmony expects genes x cells
meta_data = adata.obs[['sample']].copy()

# Run Harmony
harmony_out = hm.run_harmony(data_mat, meta_data, vars_use=['sample'], verbose=True)

# Assign the corrected PCs back to adata
adata.obsm['X_pca_harmony'] = harmony_out.Z_corr  # Already in cells x PCs format from harmonypy

# Now compute neighbors using the Harmony-corrected representation
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.tl.leiden(adata)

# Debug: check the shape of Harmony output
print('harmony_out object:', harmony_out)
print('harmony_out.Z_corr shape:', harmony_out.Z_corr.shape)
print('harmony_out.Z_corr.T shape:', harmony_out.Z_corr.T.shape)
print('Expected shape for adata.obsm:', (adata.n_obs, harmony_out.Z_corr.shape[0]))

# Remove prefix from gene names (row names)
# The prefixes are: GRCh38_ (for human genes) and GRCm39_ (for mouse genes)
adata.var_names = adata.var_names.str.replace(r'^(GRCh38_|GRCm39_)', '', regex=True)

# Verify the result
print("Gene names after removing prefix:")
print(adata.var_names[:10])
print(f"\nTotal genes: {adata.n_vars}")

# Choose resolution 0.4 which should give us close to 10 clusters
sc.tl.leiden(adata, resolution=0.1)

# Based on step 30 DEG analysis, classify clusters by species
# Human genes: ALL CAPS (e.g., EPCAM, COL1A1)
# Mouse genes: lowercase/mixed case (e.g., Picrh, Pigp3, Sparc, Sfrp2)

cluster_species = {
    '0': 'human',    # EPCAM, PIGR, etc (human markers)
    '1': 'mouse',    # Picrh, Pigp3, Sparc (mouse gene names)
    '2': 'human',    # EPCAM, KRT19 (human epithelial)
    '3': 'mouse',    # Mt-Nd4l, Mt-Nd5 (mouse mitochondrial)
    '4': 'human',    # PECAM1, VWF (human endothelial)
    '5': 'mouse',    # Cd14, Lyz (mouse immune)
    '6': 'mouse',    # Lyz, Aif1 (mouse monocytes)
    '7': 'mouse',    # Gzma, Gnly (mouse NK/T cells)
    '8': 'human',    # CD68, CD14 (human macrophages)
    '9': 'mouse',    # Acta2, Myh11 (mouse smooth muscle)
    '10': 'mouse'    # Col1a1, Postn (mouse fibroblasts)
}

# Add species column to adata
adata.obs['species'] = adata.obs['leiden'].map(cluster_species)

# Print summary
print("Species Annotation Added to adata.obs['species']")
print("\nCluster-Species Mapping:")
print("-" * 50)
for cluster, species in sorted(cluster_species.items()):
    n_cells = (adata.obs['leiden'] == cluster).sum()
    print(f"Cluster {cluster}: {species.upper():6s} ({n_cells:5d} cells)")

print("\n" + "=" * 50)
print(f"Total HUMAN cells: {(adata.obs['species'] == 'human').sum()}")
print(f"Total MOUSE cells: {(adata.obs['species'] == 'mouse').sum()}")
print("=" * 50)


# Create side-by-side UMAP plots
fig, axes = plt.subplots(1, 2, figsize=(8.5, 3))

# UMAP colored by sample
sc.pl.umap(adata, color='sample', ax=axes[0], title='UMAP - Colored by Sample', show=False)

# UMAP colored by leiden clusters (0.1 resolution)
sc.pl.umap(adata, color='layer2', ax=axes[1], title='UMAP - Colored by Leiden Clusters (Resolution 0.1)', show=False)

plt.tight_layout()
plt.show()

col = "layer2"

# Ensure categorical
adata.obs[col] = adata.obs[col].astype("category")

# (Optional safety) make sure the target label exists as a category
if "Tumor (Human)" not in adata.obs[col].cat.categories:
    adata.obs[col] = adata.obs[col].cat.add_categories(["Tumor (Human)"])

# Merge: Tumor -> Tumor (Human)
adata.obs.loc[adata.obs[col] == "Tumor", col] = "Tumor (Human)"

# Remove the unused "Tumor" category
adata.obs[col] = adata.obs[col].cat.remove_unused_categories()

col = "layer2"
s = adata.obs[col].astype("category")

cats = list(s.cat.categories)

mouse = [c for c in cats if "(Mouse)" in c]
human = [c for c in cats if "(Human)" in c]
other = [c for c in cats if ("(Mouse)" not in c and "(Human)" not in c)]

new_order = mouse + human + other

adata.obs[col] = s.cat.reorder_categories(new_order, ordered=True)

# Create dotplot with representative genes for layer1 annotation

genes = ['EPCAM','GPC3','PTPRC','CD3E','CD8A','CD4',
         'CD14','CSF1R','FCGR3A','CD68',
         'CD19','Ptprc','Cd14','C1qc','Pecam1']
# Create dotplot
sc.pl.dotplot(adata, var_names=genes, 
              groupby='layer2')
plt.tight_layout()
plt.show()

col = "layer2"

mask = adata.obs[col].isin(["Tumor (Human)"])
adata_tumor = adata[mask].copy()

import scanpy as sc

# If your subset is a view, make it a real object
adata_tumor = adata_tumor.copy()

# (Optional but usually recommended) re-run HVGs on the subset
# sc.pp.highly_variable_genes(adata_tumor, n_top_genes=2000, batch_key="sample")

# Re-run PCA on the subset
sc.tl.pca(adata_tumor)

# Recompute neighbors graph (on X_pca by default)
sc.pp.neighbors(adata_tumor, n_neighbors=15, n_pcs=30)  # tune n_pcs / n_neighbors as needed

# Recompute UMAP (optional but typical)
sc.tl.umap(adata_tumor)

sc.tl.leiden(adata_tumor, key_added="leiden_recluster", resolution=0.3, flavor="igraph", n_iterations=2)
# Create side-by-side UMAP plots
fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))

# UMAP colored by sample
sc.pl.umap(adata_tumor, color='sample', ax=axes[0], title='UMAP - Colored by Sample', show=False)

# UMAP colored by leiden clusters (0.1 resolution)
sc.pl.umap(adata_tumor, color='leiden_recluster', ax=axes[1], title='UMAP - Colored by celltypes Clusters (Resolution 0.1)', show=False)

plt.tight_layout()
plt.show()

import pandas as pd

# Make sure leiden_recluster is treated as a string key (covers "0" vs 0 cases)
lr = adata_tumor.obs["leiden_recluster"].astype(str)

mapping = {
    "2": "proliferative",
    "3": "proliferative",
    "4": "proliferative",
    "5": "proliferative",
    "6": "EMT",
    "7": "secretory",
    "0": "progenitor",
    "1": "progenitor",
    "8": "progenitor",
    "9": "progenitor",
}

adata_tumor.obs["celltype"] = lr.map(mapping).fillna("other")
adata_tumor.obs["celltype"] = pd.Categorical(
    adata_tumor.obs["celltype"],
    categories=["progenitor", "proliferative", "EMT", "secretory"],
    ordered=True
)
sc.set_figure_params(dpi=100, dpi_save=300, format="png", ipython_format="retina")

import scanpy as sc


phenotype = {
    'Oncofetal':["AFP","GPC3","LIN28B","H19","SALL4", "IGF2BP3","IGF2BP1","DLK1","CLDN6","TRIM71"],
    'Hypoxia':['HIF1A','VEGFA'],
    
               
              }


marker_dict = {
    'Progenitor':['EPCAM','CD24','SOX9'],
    'Prolif':["MKI67", "TOP2A", "CCNB1",'CDK1'],
    'EMT': ["CDH2", "FN1", "ZEB1", "SERPINE1", "SNAI1"],
    'secretory':['APOA1','APOA2','ALB'],
               
                  
}

sc.pl.dotplot(
    adata_tumor,
    marker_dict,
    groupby="celltype",
    standard_scale="var",
)

sc.pl.dotplot(
    adata_tumor,
    bad,
    groupby="celltype",
    standard_scale="var",
)

import scanpy as sc

sc.pl.matrixplot(
    adata_tumor,
    bad,                   # list (or dict) of genes
    groupby="celltype",
    standard_scale="var",  # same idea as in dotplot
    cmap="viridis",
    dendrogram=False,       # optional
)

import scanpy as sc

import scanpy as sc

sc.pl.dotplot(
    adata_tumor,
    marker_dict,
    groupby="celltype",
    standard_scale="var",
    swap_axes=True,   # <- make it “vertical” [web:380]
)


sc.pl.matrixplot(
    adata_tumor,
    bad,
    groupby="celltype",
    standard_scale="var",
    cmap="viridis",
    dendrogram=False,
    swap_axes=True,   # <- makes groupby categories on x-axis [web:228]
)

import pandas as pd
import matplotlib.pyplot as plt

sample_col   = "sample"
celltype_col = "celltype"

# rows = sample, cols = celltype, values = % within each sample
tab_pct = (
    pd.crosstab(
        adata_tumor.obs[sample_col],
        adata_tumor.obs[celltype_col],
        normalize="index",
        dropna=False,
    )
    .mul(100)
)

# Optional: enforce stacking order if celltype is ordered categorical
if pd.api.types.is_categorical_dtype(adata_tumor.obs[celltype_col]):
    tab_pct = tab_pct.reindex(columns=adata_tumor.obs[celltype_col].cat.categories)

fig, ax = plt.subplots(figsize=(0.7 * tab_pct.shape[0] + 3, 5))

tab_pct.plot(
    kind="bar",          # vertical stacked bars [web:351]
    stacked=True,
    ax=ax,
    edgecolor="none",
    linewidth=0,
)

ax.set_xlabel(sample_col)
ax.set_ylabel("Percent of cells")
ax.set_ylim(0, 100)

ax.grid(False)
leg = ax.get_legend()
handles = leg.legend_handles
labels  = [t.get_text() for t in leg.get_texts()]
leg.remove()

ax.legend(
    handles[::-1], labels[::-1],
    title=celltype_col, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False
)

plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.show()

import scanpy as sc

oncofetal_genes = [
    "AFP","GPC3","LIN28B","H19","SALL4",
    "IGF2BP3","IGF2BP1","DLK1","CLDN6","TRIM71"
]

# Pick where genes should be looked up (adata.var_names vs adata.raw.var_names)
use_raw = adata.raw is not None  # score_genes defaults to using .raw if present when use_raw=None [web:1]
varnames = adata.raw.var_names if use_raw else adata.var_names

genes_present = [g for g in oncofetal_genes if g in varnames]
genes_missing = [g for g in oncofetal_genes if g not in varnames]

# Compute per-cell module score (written to adata.obs["oncofetal_score"]) [web:1]
sc.tl.score_genes(
    adata,
    gene_list=genes_present,
    score_name="oncofetal_score",
    use_raw=use_raw,     # set explicitly to avoid surprises [web:1]
    # layer="log1p",      # optionally score from a specific layer instead of .X/.raw [web:1]
    # ctrl_size=50,       # default reference genes per bin [web:1]
    # n_bins=25,          # default expression bins [web:1]
)

print("Missing genes:", genes_missing)
print(adata.obs[["oncofetal_score"]].head())

import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# Compute p (from your data)
g1 = adata.obs[adata.obs['sample']=='carm']['oncofetal_score']
g2 = adata.obs[adata.obs['sample']=='carmdr']['oncofetal_score']
_, p = mannwhitneyu(g1, g2)

def sig_star(p):
    if p < 1e-4: return '****'
    elif p < 1e-3: return '***'
    elif p < 1e-2: return '**'
    elif p < 0.05: return '*'
    else: return 'ns'

star = sig_star(p)
print(f"carm vs carmdr: p={p}, significance: {star}")

# Your plot
ax = sc.pl.violin(adata, 'oncofetal_score', groupby='sample', rotation=45, stripplot=False, show=False)
ax.grid(False)

# Add star manually (adjust y_pos to max + offset)
y_max = max(g1.max(), g2.max())
ax.text(0.5, y_max * 1.05, star, ha='center', va='bottom', fontsize=16, fontweight='bold')
ax.set_ylim(top=y_max * 1.15)
plt.show()
