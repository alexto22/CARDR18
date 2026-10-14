"""Explicit reconstruction, NOT recovery of historical preprocessing.

Chunked mixed-human/mouse raw UMI -> per-cell total 10,000 -> natural log1p.
No cell filtering or inferred annotations. Scores follow the expression-bin
algorithm documented in Scanpy 1.9.8 tools/_score_genes.py, implemented with
installed numpy/pandas/h5py because Scanpy/AnnData are unavailable locally.
Historical Scanpy version, normalization, raw layer and cluster-7 membership
are unknown. The resulting HDF5 is a documented COO store, NOT an .h5ad.
"""
import argparse
import gzip
import hashlib
import itertools
import json
import platform
import re
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SIGNATURE = ['AFP', 'GPC3', 'LIN28B', 'H19', 'SALL4', 'IGF2BP3',
             'IGF2BP1', 'DLK1', 'CLDN6', 'TRIM71']
CHUNK = 250000


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def unique_names(names):
    occupied = set(names)
    seen = set()
    result = []
    counters = {}
    for name in names:
        if name not in seen:
            result.append(name)
            seen.add(name)
            continue
        number = counters.get(name, 0)
        while True:
            number += 1
            candidate = name + '-' + str(number)
            if candidate not in occupied:
                break
        counters[name] = number
        occupied.add(candidate)
        result.append(candidate)
    return result


def text_gzip(path):
    with gzip.open(path, 'rt') as f:
        return [s.rstrip('\n') for s in f]


def matrix_chunks(path):
    # Reading to EOF checks gzip CRC/length. Never allocate a dense matrix.
    with gzip.open(path, 'rt') as f:
        for line in f:
            if not line.startswith('%'):
                dims = tuple(map(int, line.split()))
                break
        yield dims
        while True:
            lines = list(itertools.islice(f, CHUNK))
            if not lines:
                break
            a = np.fromstring(''.join(lines), dtype=np.int64, sep=' ')
            if a.size != 3 * len(lines):
                raise ValueError('Invalid Matrix Market entries')
            yield a.reshape(-1, 3)


def initialize_group(store, name, barcodes, genes, nnz):
    g = store.create_group(name)
    dt = h5py.string_dtype('utf-8')
    g.create_dataset('barcodes', data=np.array(barcodes, dtype=object), dtype=dt)
    g.create_dataset('gene_ids', data=np.array(genes, dtype=object), dtype=dt)
    for key in ['row_gene_index', 'column_cell_index', 'raw_umi']:
        g.create_dataset(key, shape=(nnz,), dtype='uint32',
                         chunks=(min(CHUNK, nnz),), compression='gzip',
                         compression_opts=1)
    g.attrs['orientation'] = 'COO records: gene rows, cell columns, zero-based'
    return g


def write_entries(g, start, rows, cols, values, totals, n_genes):
    if (np.any(rows < 0) or np.any(rows >= n_genes)
            or np.any(cols < 0) or np.any(cols >= len(totals))
            or np.any(values < 0) or np.any(values > np.iinfo('uint32').max)):
        raise ValueError('Invalid raw count/index')
    end = start + len(values)
    g['row_gene_index'][start:end] = rows
    g['column_cell_index'][start:end] = cols
    g['raw_umi'][start:end] = values
    totals += np.bincount(cols, weights=values, minlength=len(totals))


def write_provenance(out, hashes, car_shape, dr_shape, car_barcodes, dr_barcodes, controls):
    evidence = {'status': 'reconstruction; not exact historical reproduction',
        'input_files': hashes, 'python': platform.python_version(),
        'versions': {'numpy': np.__version__, 'pandas': pd.__version__, 'h5py': h5py.__version__},
        'retained_cells': {'carm': len(car_barcodes), 'carmdr': len(dr_barcodes)},
        'historical_plot_cells': {'carm': 25094, 'carmdr': 18704},
        'raw_dimensions': {'carm': car_shape, 'carmdr': dr_shape},
        'signature': SIGNATURE, 'control_genes': sorted(controls),
        'scoring': {'algorithm_reference': 'Scanpy 1.9.8 tools/_score_genes.py',
            'algorithm_url': 'https://raw.githubusercontent.com/scverse/scanpy/1.9.8/scanpy/tools/_score_genes.py',
            'n_bins': 25, 'ctrl_size': 50, 'random_state': 0,
            'gene_pool': 'all 157230 measured mixed-species Gene Expression features',
            'expression': 'explicitly reconstructed log1p matrix; no historical .raw recovered'},
        'assumptions': ['Normalize each cell across all measured genes to 10000 UMI, then natural log1p',
            'Retain all raw cells: historical cluster-7 filtering unavailable',
            'No Harmony/clustering/annotation needed or rerun for the score',
            'No human-only or tumour-only selection; annotations unavailable'],
        'limitations': ['Historical preprocessing and Scanpy version unknown',
            'Original plot pooled mixed cell types; manuscript says cancer cells',
            'Animal/biological replicate metadata absent; cell-level p is descriptive'],
        'integrity': 'CAR-M gzip stream read to EOF/CRC; MTX dimensions and entry count checked; '
                     'H5 CSC pointers, feature alignment, dimensions and nonnegative counts checked'}
    (out / 'reconstruction_provenance.json').write_text(json.dumps(evidence, indent=2, default=int))
    print('Completed', out.resolve(), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path,
                        default=Path('hcc_raw_reconstruction'))
    parser.add_argument('--resume-completed', action='store_true',
                        help='Write missing metadata for a completed store, without recomputing')
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    store_path = out / 'processed_mixed_cells_coo.h5'
    if store_path.exists() and args.resume_completed:
        with h5py.File(store_path, 'r') as f:
            if f.attrs['status'] != 'complete':
                raise ValueError('Cannot resume an incomplete store')
            car_barcodes = list(f['carm/barcodes'][:])
            dr_barcodes = list(f['carmdr/barcodes'][:])
            n_genes = len(f['scoring_gene_names'])
            car_shape = (n_genes, len(car_barcodes), len(f['carm/raw_umi']))
            dr_shape = (n_genes, len(dr_barcodes))
        hashes = json.loads((out / 'input_hashes.json').read_text())
        audit = pd.read_csv(out / 'scoring_gene_audit.csv')
        controls = set(audit.loc[audit['is_control'], 'gene'])
        write_provenance(out, hashes, car_shape, dr_shape, car_barcodes, dr_barcodes, controls)
        return
    if store_path.exists():
        raise FileExistsError('Preserving existing store; choose a new output directory')
    inputs = [args.input_root / 'carm' / name for name in
              ['matrix.mtx.gz', 'barcodes.tsv.gz', 'features.tsv.gz']]
    inputs.append(args.input_root / 'cardrm/filtered_feature_bc_matrix.h5')
    print('Hashing source files', flush=True)
    hashes = {str(p.relative_to(args.input_root)): {'bytes': p.stat().st_size,
                               'sha256': sha256(p)} for p in inputs}
    (out / 'input_hashes.json').write_text(json.dumps(hashes, indent=2))
    car_barcodes = text_gzip(inputs[1])
    feature_rows = [s.split('\t') for s in text_gzip(inputs[2])]
    gene_ids = [s[0] for s in feature_rows]
    symbols = [s[1] for s in feature_rows]
    if any(len(s) != 3 or s[2] != 'Gene Expression' for s in feature_rows):
        raise ValueError('Non-Gene-Expression features require explicit filtering')
    with h5py.File(inputs[3], 'r') as raw:
        m = raw['matrix']
        dr_shape = tuple(map(int, m['shape'][:]))
        dr_ids = [x.decode() for x in m['features/id'][:]]
        dr_symbols = [x.decode() for x in m['features/name'][:]]
        dr_barcodes = [x.decode() for x in m['barcodes'][:]]
        if dr_ids != gene_ids or dr_symbols != symbols:
            raise ValueError('Gene order differs: explicit union alignment required')
        if any(x.decode() != 'Gene Expression'
               for x in m['features/feature_type'][:]):
            raise ValueError('Unexpected H5 feature type')
        dr_ptr = m['indptr'][:]
        dr_nnz = len(m['data'])
        assert dr_shape == (len(gene_ids), len(dr_barcodes))
        assert len(dr_ptr) == len(dr_barcodes) + 1
        assert dr_ptr[0] == 0 and dr_ptr[-1] == dr_nnz
        assert np.all(np.diff(dr_ptr) >= 0)
        assert len(m['indices']) == dr_nnz
        # Same name handling sequence as the recorded notebook.
        names = unique_names([re.sub(r'^(GRCh38_|GRCm39_)', '', x)
                              for x in unique_names(symbols)])
        pd.DataFrame({'gene_index': np.arange(len(names)), 'gene_id': gene_ids,
                      'original_symbol': symbols, 'scoring_name': names}).to_csv(
                          out / 'gene_metadata.csv', index=False)
        chunks = matrix_chunks(inputs[0])
        car_shape = next(chunks)
        assert car_shape[:2] == (len(gene_ids), len(car_barcodes))
        n_cells = len(car_barcodes) + len(dr_barcodes)
        gene_sum = np.zeros(len(names), dtype=np.float64)
        totals = {}
        with h5py.File(store_path, 'w') as store:
            store.attrs['status'] = 'incomplete'
            store.attrs['reconstruction'] = True
            store.attrs['normalization'] = 'total RNA UMI to 10000, then natural log1p'
            store.attrs['filtering'] = 'none; historical cluster-7 assignment unavailable'
            dt = h5py.string_dtype('utf-8')
            store.create_dataset('scoring_gene_names',
                                 data=np.array(names, dtype=object), dtype=dt)
            car = initialize_group(store, 'carm', car_barcodes, gene_ids, car_shape[2])
            totals['carm'] = np.zeros(len(car_barcodes), dtype=np.float64)
            start = 0
            print('Validating/streaming CAR-M matrix', car_shape, flush=True)
            for a in chunks:
                write_entries(car, start, a[:, 0]-1, a[:, 1]-1, a[:, 2],
                              totals['carm'], len(names))
                start += len(a)
                if start % 5000000 == 0:
                    print('CAR-M entries', start, flush=True)
            assert start == car_shape[2], 'Truncated or excess MTX records'
            dr = initialize_group(store, 'carmdr', dr_barcodes, gene_ids, dr_nnz)
            totals['carmdr'] = np.zeros(len(dr_barcodes), dtype=np.float64)
            print('Validating/streaming DR18 H5', dr_shape, flush=True)
            for start in range(0, dr_nnz, CHUNK):
                end = min(start + CHUNK, dr_nnz)
                columns = np.searchsorted(dr_ptr, np.arange(start, end),
                                          side='right') - 1
                write_entries(dr, start, m['indices'][start:end], columns,
                              m['data'][start:end], totals['carmdr'], len(names))
            print('Normalization and gene means', flush=True)
            for label in ['carm', 'carmdr']:
                g = store[label]
                if np.any(totals[label] <= 0):
                    raise ValueError('Zero-UMI cell: no silent filtering permitted')
                g.create_dataset('total_raw_umi', data=totals[label])
                nnz = len(g['raw_umi'])
                norm = g.create_dataset('log1p_total_10000', shape=(nnz,),
                    dtype='float64', chunks=(min(CHUNK, nnz),),
                    compression='gzip', compression_opts=1)
                for start in range(0, nnz, CHUNK):
                    end = min(start + CHUNK, nnz)
                    rows = g['row_gene_index'][start:end]
                    cols = g['column_cell_index'][start:end]
                    values = np.log1p(g['raw_umi'][start:end].astype(np.float64)
                                     * 10000 / totals[label][cols])
                    norm[start:end] = values
                    gene_sum += np.bincount(rows, weights=values,
                                            minlength=len(names))
            avg = pd.Series(gene_sum / n_cells, index=names)
            bins = avg.rank(method='min') // int(np.round(len(avg) / 24))
            missing = [g for g in SIGNATURE if g not in names]
            if missing:
                raise ValueError('Missing signature genes: ' + str(missing))
            rng = np.random.RandomState(0)
            controls = set()
            for cut in np.unique(bins.loc[SIGNATURE]):
                candidates = np.array(bins[bins == cut].index)
                rng.shuffle(candidates)
                controls.update(candidates[:50])
            controls -= set(SIGNATURE)
            control_indices = [i for i, g in enumerate(names) if g in controls]
            weights = np.zeros(len(names), dtype=np.float64)
            weights[control_indices] = -1 / len(controls)
            signature_indices = [names.index(g) for g in SIGNATURE]
            weights[signature_indices] = 1 / len(SIGNATURE)
            pd.DataFrame({'gene': names, 'mean_log1p': avg.values,
                          'expression_bin': bins.values,
                          'is_signature': [g in SIGNATURE for g in names],
                          'is_control': [g in controls for g in names]}).to_csv(
                              out / 'scoring_gene_audit.csv', index=False)
            frames = []
            print('Calculating scores with', len(controls), 'control genes', flush=True)
            for label, barcodes in [('carm', car_barcodes), ('carmdr', dr_barcodes)]:
                g = store[label]
                score = np.zeros(len(barcodes), dtype=np.float64)
                for start in range(0, len(g['raw_umi']), CHUNK):
                    end = min(start + CHUNK, len(g['raw_umi']))
                    rows = g['row_gene_index'][start:end]
                    cols = g['column_cell_index'][start:end]
                    values = g['log1p_total_10000'][start:end]
                    score += np.bincount(cols, weights=values * weights[rows],
                                         minlength=len(barcodes))
                g.create_dataset('oncofetal_score', data=score)
                frames.append(pd.DataFrame({'cell_id': [label + ':' + b for b in barcodes],
                    'barcode': barcodes, 'sample': label,
                    'total_raw_umi': totals[label], 'oncofetal_score': score}))
            pd.concat(frames, ignore_index=True).to_csv(out / 'cell_scores.csv', index=False)
            store.attrs['status'] = 'complete'
    write_provenance(out, hashes, car_shape, dr_shape, car_barcodes, dr_barcodes, controls)


if __name__ == '__main__':
    main()
