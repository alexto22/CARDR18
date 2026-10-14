# HCC mixed-cell oncofetal reconstruction

This is a **new, explicitly specified reconstruction**, not recovery of the historical processed object or an exact reproduction of the historical figure. It is separate from the existing manuscript notebooks in this repository.

## Processed object and source data

[Restricted Dropbox reconstruction folder](https://www.dropbox.com/scl/fo/36j1v2sp00q480jam5ph8/h?rlkey=vqdpyql5bf1j7yky2zegmgiyv&dl=0)

Download `processed_mixed_cells_coo.h5` from that folder: 578,597,688 bytes (551.79 MiB). This is a custom sparse HDF5 COO store, **not AnnData / `.h5ad`**. The folder also contains `cell_scores.csv`, `reconstruction_provenance.json`, gene metadata, scoring audits and validation reports. Access is restricted to people already granted Dropbox access; posting this link does not grant public access. Contact the repository owner if access is needed. No object binary or raw matrices are committed here.

## What was reconstructed

- CAR-M: 25,296 cells; CAR-M_DR18: 20,026 cells; all 157,230 shared mixed-human/mouse Gene Expression features retained.
- All cells present in the supplied matrices are retained; there is no additional cell-type or species selection.
- Per-cell total normalization to 10,000 UMI, then natural-log `log1p`.
- Signature: AFP, GPC3, LIN28B, H19, SALL4, IGF2BP3, IGF2BP1, DLK1, CLDN6, TRIM71.
- Local implementation of the [Scanpy 1.9.8 score_genes algorithm](https://github.com/scverse/scanpy/blob/1.9.8/scanpy/tools/_score_genes.py): 25 expression bins, 50 candidate control genes per occupied bin, random seed 0; all measured features form the control pool. The original analysis's Scanpy version is unknown.
- CAR-M is red (`#D62728`); CAR-M_DR18 is blue (`#1F77B4`).

Historical normalization, raw layer and cluster-7 membership were not recovered. The historical plotted populations were 25,094 and 18,704 cells, respectively. This reconstruction does not implement that historical filtering and must not be described as tumour-only evidence. Biological replicate metadata is unavailable: the cell-level Mann–Whitney test is descriptive and does not establish animal-level replication. A reported floating-point p-value of zero means numerical underflow, not a mathematically zero probability.

## Plot only

Download `cell_scores.csv` and `reconstruction_provenance.json` into the same local folder. This path does not load matrices, normalize or rescore cells.

```sh
python regenerate_oncofetal_plot.py /path/to/cell_scores.csv --output-dir ./plots
```

Outputs: PDF, SVG and 600-dpi PNG, score summary, statistics JSON and caption. Input cell counts, sample labels, finite scores and unique cell IDs are checked.

## Reconstruct from input matrices

Prepare this layout (the spelling `cardrm` matches the source layout):

```text
HCC/
  carm/matrix.mtx.gz
  carm/barcodes.tsv.gz
  carm/features.tsv.gz
  cardrm/filtered_feature_bc_matrix.h5
```

```sh
python reconstruct_hcc_oncofetal_from_raw.py --input-root /path/to/HCC --output-dir ./reconstructed
python regenerate_oncofetal_plot.py ./reconstructed/cell_scores.csv
```

Use a new output directory: the script refuses to overwrite an existing store. Do not run Python with `-O`, which disables input assertions. `--resume-completed` only regenerates missing provenance for a completed store; it does not resume interrupted matrix processing. The public script requires an explicit input root and records relative input paths, removing the original machine-specific default; the numerical method is unchanged.

## Read the processed object

Each `carm` / `carmdr` group contains `barcodes`, `gene_ids`, zero-based `row_gene_index` and `column_cell_index`, `raw_umi`, `log1p_total_10000`, `total_raw_umi` and `oncofetal_score`. Top-level `scoring_gene_names` supplies the scoring feature names. Matrix orientation is genes × cells. Use chunked access for the large coordinate arrays; do not densify the full matrix.

```python
import h5py
with h5py.File('processed_mixed_cells_coo.h5', 'r') as f:
    assert f.attrs['status'] == 'complete'
    for label in ('carm', 'carmdr'):
        g = f[label]
        scores = g['oncofetal_score'][:]
        print(label, len(g['barcodes']), len(g['gene_ids']), scores.mean())
        # Example sparse coordinate slice, not a dense matrix:
        rows = g['row_gene_index'][:1000]
        cols = g['column_cell_index'][:1000]
        values = g['log1p_total_10000'][:1000]
```

## Dependencies and verification

Dependencies: NumPy, pandas, h5py, SciPy, Matplotlib and seaborn. Scanpy and AnnData are not required by these scripts. The original reconstruction ran with Python 3.8.5, NumPy 1.21.2, pandas 1.3.5, h5py 2.10.0; original plotting additionally used SciPy 1.7.3, Matplotlib 3.5.0, seaborn 0.11.2. Use those versions in an isolated compatible environment when matching the original reconstruction's rendering.

Publication checks: both scripts compile; the plot-only script ran on the supplied score CSV with Python 3.12.14, NumPy 2.3.5, pandas 2.2.3, SciPy 1.17.0, Matplotlib 3.10.8 and seaborn 0.13.2. It reproduced cell counts and Mann–Whitney U = 480,293,138.5; newer seaborn emitted deprecation warnings. Full raw reconstruction was not rerun during publication; this publication environment lacks h5py, so the raw-reconstruction CLI was not runtime-tested here. Historical preprocessing and exact historical reproduction remain unverified.

Score CSV SHA-256: `50100d34391cf059be60de307b5652080e7b24e3a636296069cb4d445e89deb5`.
