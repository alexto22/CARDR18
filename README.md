# IL-18BP-resistant CAR-macrophages for solid tumor therapy

CAR-macrophages (CAR-Ms) are a promising solid tumor therapy, but their activity is limited by immunosuppressive feedback within the tumor microenvironment. To address this limitation, we defined the immune checkpoint landscape of the tumor microenvironment and found upregulation of IL-18BP in CAR-Ms. Genetic deletion of IL-18BP confirmed its detrimental role in CAR-Ms. To overcome IL-18BP-mediated IL-18 neutralization, we engineered CAR-Ms to express decoy-resistant IL-18 (DR-18). In syngeneic immunocompetent models, DR-18-armored CAR-M improved tumor control and long-term survival compared with conventional CAR-M and IL-18-armored CAR-M. Mechanistically, DR-18 armoring retained the M1-like state of intratumoral CAR-M and activated NK- and T-cells, achieving the reduction of oncofetal signatures and tumor-initiating cells. Finally, we found that DR-18 applies to CAR-NK. Together, these findings define IL-18BP as a barrier to CAR-M therapy, which can be overcome by DR-18.

## Download links for processed files

- [combined_macrophages_normalized.rds](https://www.dropbox.com/scl/fi/e9oufws4j6cu9ojqo7h0r/combined_macrophages_normalized.rds?rlkey=yj81lu5aa5a8ihe1yp83vnozo&dl=1)
- [three_sample_macrophage_CAR_renamed.rds](https://www.dropbox.com/scl/fi/dvwpl31i598gtqsgxbdpn/three_sample_macrophage_CAR_renamed.rds?rlkey=t1sqluhbmemf3nm8cv9i4lfu3&dl=1)

## Analysis scripts in this repository

```text
CARDR18/
├── Mac_v3.ipynb             # Figure 1f–i
├── CARM_DR18_v2.ipynb        # Figure 4f–h
├── my_list_barplots.R
└── signature_barplots.R
```

Updated: 12 October 2026

## HCC analysis

- [Processed HCC files](https://www.dropbox.com/scl/fo/36j1v2sp00q480jam5ph8/h?rlkey=vqdpyql5bf1j7yky2zegmgiyv&dl=0): `processed_mixed_cells_coo.h5` and `cell_scores.csv`. Dropbox access is restricted to existing authorized viewers. The HDF5 file is a sparse COO store, not an AnnData `.h5ad` file.
- [HCC_analysis.py](HCC_analysis.py): prepare the mixed-cell matrix and calculate oncofetal scores.
- [HCC_plot.py](HCC_plot.py): plot the saved scores, with CAR-M in red and CAR-M_DR18 in blue. Keep the accompanying provenance JSON from Dropbox beside `cell_scores.csv`.

The analysis retains all supplied cells (25,296 CAR-M; 20,026 CAR-M_DR18) and mixed-human/mouse features, normalizes each cell to 10,000 total UMI and applies natural-log1p. Scores use the Scanpy 1.9.8 expression-bin method. Historical cluster-7 filtering is not applied; these are mixed-cell results, with descriptive cell-level statistics rather than biological-replicate inference.

```sh
python HCC_analysis.py --input-root /path/to/HCC --output-dir ./HCC_analysis_output
python HCC_plot.py /path/to/cell_scores.csv --output-dir ./plots
```
