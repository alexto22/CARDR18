# IL-18BP-resistant CAR-macrophages for solid tumor therapy

CAR-macrophages (CAR-Ms) are a promising solid tumour therapy, but their activity is limited by immunosuppressive feedback within the tumour microenvironment. This study identifies IL-18 binding protein (IL-18BP) as a barrier to CAR-M activity and engineers CAR-Ms to express decoy-resistant IL-18 (DR-18). DR-18 armouring improves tumour control, sustains inflammatory macrophage programmes and activates NK- and T-cell responses in preclinical models.

This repository contains single-cell RNA-sequencing analysis code and links to processed macrophage objects accompanying the manuscript.

## Processed data

Download the processed Seurat objects below and keep them alongside the notebooks.

- [combined_macrophages_normalized.rds](https://www.dropbox.com/scl/fi/e9oufws4j6cu9ojqo7h0r/combined_macrophages_normalized.rds?rlkey=yj81lu5aa5a8ihe1yp83vnozo&dl=1) — normalised Seurat object for the Figure 1 comparison of CAR-M, tumour-associated macrophages and adjacent-normal macrophages; 6,082 cells.
- [three_sample_macrophage_CAR_renamed.rds](https://www.dropbox.com/scl/fi/dvwpl31i598gtqsgxbdpn/three_sample_macrophage_CAR_renamed.rds?rlkey=t1sqluhbmemf3nm8cv9i4lfu3&dl=1) — Seurat object used for the Figure 4 CAR-M versus DR-18 CAR-M comparison; 438 cells. Despite its filename, the comparison uses the `CAR` and `CAR_DR18` labels.

## Analysis scripts in this repository

```text
CARDR18/
├── Fig1_IL18BP.r                 # Existing IL-18BP macrophage comparison
├── Fig3_CARM.r                   # Existing CAR-M characterisation plots
├── Fig4_hCC.py                   # Existing hCC single-cell analysis
├── Mac_v3.ipynb                  # Figure 1f–i: CAR-M, tumour and adjacent-normal macrophages
├── CARM_DR18_v2.ipynb             # Figure 4f–h: CAR-M versus DR-18 CAR-M
├── my_list_barplots.R            # Gene-panel plotting helper for CARM_DR18_v2.ipynb
└── signature_barplots.R          # Signature-panel plotting helper for CARM_DR18_v2.ipynb
```

## Running the macrophage notebooks

Use Jupyter with an R kernel (IRkernel), Seurat and SeuratObject 5 or later, ggplot2, patchwork, Matrix and IRdisplay. Place the processed RDS files alongside the notebooks and helper scripts, then run each notebook from that directory in cell order.

`Mac_v3.ipynb` loads `combined_macrophages_normalized.rds` by default (`force_rebuild <- FALSE`). Its optional rebuild path additionally requires the original HCC reference object `HCCF1F4.RData`, which is not included here. `CARM_DR18_v2.ipynb` loads `three_sample_macrophage_CAR_renamed.rds` and sources the two plotting helpers.

The CAR-M/DR-18 notebook writes its tables and provenance records to `results/carm_dr18_two_panels`. Its exploratory tests treat pooled cells as observations and do not establish biological-replicate significance. Consult the notebook for normalisation, gene panels and statistical details.

The distributed notebooks have cleared cached outputs and portable output paths. The three existing R/Python scripts are retained unchanged and may require additional original inputs and path configuration. Their filenames are retained from earlier analyses. Full end-to-end execution has not been independently validated for this release.

Updated: 12 October 2026
