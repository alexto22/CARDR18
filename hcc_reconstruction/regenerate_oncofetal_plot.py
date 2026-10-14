"""Plot only reconstructed HCC mixed-cell scores; no raw loading or rescoring.
Input: cell_scores.csv plus reconstruction_provenance.json from the separate
reconstruct_hcc_oncofetal_from_raw.py. All outputs are labelled reconstruction.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scores', type=Path)
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scipy
    import seaborn as sns
    from scipy.stats import mannwhitneyu

    out = args.output_dir or args.scores.parent
    out.mkdir(parents=True, exist_ok=True)
    provenance = json.loads((args.scores.parent / 'reconstruction_provenance.json').read_text())
    obs = pd.read_csv(args.scores)
    if set(obs['sample']) != {'carm', 'carmdr'}:
        raise ValueError('Expected sample labels carm/carmdr')
    if not np.isfinite(obs['oncofetal_score'].to_numpy()).all():
        raise ValueError('Non-finite scores')
    counts = obs['sample'].value_counts().to_dict()
    if counts != provenance['retained_cells']:
        raise ValueError('Score population differs from reconstruction provenance')
    if obs['cell_id'].duplicated().any():
        raise ValueError('Duplicated cell identifier')
    obs['sample'] = pd.Categorical(obs['sample'], categories=['carm', 'carmdr'])
    first = obs.loc[obs['sample'] == 'carm', 'oncofetal_score']
    second = obs.loc[obs['sample'] == 'carmdr', 'oncofetal_score']
    statistic, p = mannwhitneyu(first, second, alternative='two-sided',
                               use_continuity=True, method='auto')
    star = next((s for t, s in [(1e-4, '****'), (1e-3, '***'),
                               (1e-2, '**'), (0.05, '*')] if p < t), 'ns')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
        'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'axes.linewidth': 0.75, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'svg.fonttype': 'none', 'figure.facecolor': 'white'})
    figure, ax = plt.subplots(figsize=(85/25.4, 85/25.4))
    figure.subplots_adjust(left=0.22, right=0.97, bottom=0.24, top=0.86)
    sns.violinplot(data=obs, x='sample', y='oncofetal_score',
        order=['carm', 'carmdr'], palette={'carm': '#D62728', 'carmdr': '#1F77B4'},
        saturation=1, inner=None, cut=0, scale='width', bw='scott',
        gridsize=100, linewidth=0.85, ax=ax)
    y_min = float(obs['oncofetal_score'].min())
    y_max = float(obs['oncofetal_score'].max())
    span = max(y_max-y_min, 0.1)
    ax.text(0.5, y_max + span*0.065, star, ha='center', va='bottom', fontsize=11)
    ax.set_ylim(y_min-span*0.04, y_max+span*0.20)
    ax.set_ylabel('Oncofetal Score')
    ax.set_xlabel('')
    ax.set_xticklabels(['CAR-M\n(n={:,} cells)'.format(counts['carm']),
                        'CAR-M_DR18\n(n={:,} cells)'.format(counts['carmdr'])])
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    figure.suptitle('Mixed-cell reconstruction', fontsize=9, y=0.96)
    figure.text(0.22, 0.055, 'Cell-level Mann–Whitney test\nBiological replicates unavailable',
                fontsize=8, ha='left', va='bottom')
    for extension in ['pdf', 'svg', 'png']:
        figure.savefig(out / ('oncofetal_score_reconstruction.'+extension), dpi=600)
    summary = obs.groupby('sample', observed=True)['oncofetal_score'].agg(
        ['count', 'mean', 'median', 'std', 'min', 'max'])
    summary.to_csv(out / 'score_summary.csv')
    result = {'status': provenance['status'], 'input': str(args.scores.resolve()),
        'score_csv_sha256': hashlib.sha256(args.scores.read_bytes()).hexdigest(),
        'n_cells': counts, 'mannwhitney_u': float(statistic), 'p_numeric': float(p),
        'p_display': 'numerical underflow (not mathematically zero)' if p == 0 else repr(float(p)),
        'stars': star, 'test': 'two-sided Mann-Whitney U; continuity correction; '
              'auto method (asymptotic for these large groups); tie correction; no multiplicity adjustment',
        'versions': {'python': platform.python_version(), 'numpy': np.__version__,
            'pandas': pd.__version__, 'scipy': scipy.__version__,
            'matplotlib': matplotlib.__version__, 'seaborn': sns.__version__},
        'colors': {'CAR-M': '#D62728', 'CAR-M_DR18': '#1F77B4'},
        'violin': {'bw': 'scott', 'cut': 0, 'scale': 'width', 'gridsize': 100,
            'inner': None, 'figure_mm': [85,85], 'png_dpi': 600},
        'historical_saved_p': 0.0,
        'limitation': 'Cell-level inference only. Original preprocessing/filtering unrecovered.'}
    (out / 'plot_statistics.json').write_text(json.dumps(result, indent=2))
    (out / 'figure_caption.txt').write_text(
        'Reconstructed mixed-cell oncofetal scores from two 10x libraries, retaining all '
        'raw cells (CAR-M n={:,}; CAR-M_DR18 n={:,}). Each cell was normalized to 10,000 '
        'total UMI across the shared mixed-human/mouse feature set and natural-log1p '
        'transformed. Scores subtract expression-matched control-gene means from the '
        'mean of the ten specified signature genes (25 bins, 50 candidate controls per '
        'occupied bin, random seed 0; Scanpy 1.9.8 algorithm implemented locally). '
        'Violins use Scott KDE bandwidth and equal maximum widths. Significance is '
        'the descriptive, two-sided cell-level Mann–Whitney U test; biological replicate '
        'metadata is unavailable. Historical cluster-7 removal and preprocessing could '
        'not be recovered; this is not an exact reproduction or tumour-only analysis.\n'.format(
            counts['carm'], counts['carmdr']))
    print(summary.to_string())
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
