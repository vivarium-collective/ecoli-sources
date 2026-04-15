"""
compare_datasets.py

Scatter plot and concordance report comparing TPM values between two datasets,
colored by gene category (essential × expression-adjusted).

Usage:
    uv run python analysis/compare_datasets.py
    uv run python analysis/compare_datasets.py --x vecoli_m9_glucose_minus_aas \
        --y gbw_vegas_wt_m9glc_34h_no_ssrA --n_annotate 20
    uv run python analysis/compare_datasets.py \
        --y gbw_vegas_wt_m9glc_34h_no_ssrA --report
    uv run python analysis/compare_datasets.py \
        --y gbw_vegas_wt_m9glc_34h_no_ssrA --report --log2fc_threshold -2
"""

import argparse
from datetime import date
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from pathlib import Path
from scipy import stats

HERE = Path(__file__).parent
DATA_DIR = HERE / '../data_formatted'
RESULTS_DIR = HERE / 'model_results'   # vEcoli modeling outputs (interface with vEcoli repo)
REPORTS_DIR = HERE / 'results'          # concordance reports and other analysis outputs
FIGURES_DIR = HERE / 'figures'
FIGURES_DIR.mkdir(exist_ok=True)
REPORTS_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Gene categories (essential × expression-adjusted)
# ---------------------------------------------------------------------------

PALETTE = {
    'other':                  '#cccccc',
    'essential':              '#4c78a8',
    'expression adjusted':    '#e45756',
    'essential + adjusted':   '#f58518',
}
ZORDER = {cat: i for i, cat in enumerate(PALETTE)}


def assign_category(row):
    ess = bool(row.get('is_essential', False))
    adj = bool(row.get('has_expression_adjustment', False))
    if ess and adj:
        return 'essential + adjusted'
    if ess:
        return 'essential'
    if adj:
        return 'expression adjusted'
    return 'other'


# ---------------------------------------------------------------------------
# Plotting function
# ---------------------------------------------------------------------------

def tpm_comparison_scatter(
    ax,
    x, y, labels, categories,
    xlabel, ylabel, title,
    n_annotate=15,
):
    """
    Log-log TPM comparison scatter.

    Genes with TPM=0 in one dataset are shown as off-axis triangles:
      - x=0, y>0 : left-pointing triangle at the left clip boundary
      - x>0, y=0 : downward triangle at the bottom clip boundary

    Genes with both=0 are omitted.

    Outliers (largest |log2 fold-change from identity|) among plotted genes
    are annotated with their gene symbol.
    """
    # Reset indices so boolean indexing stays consistent throughout
    x = x.reset_index(drop=True)
    y = y.reset_index(drop=True)
    labels = labels.reset_index(drop=True)
    categories = categories.reset_index(drop=True)

    both_pos  = (x > 0) & (y > 0)
    x_zero    = (x == 0) & (y > 0)   # not in x dataset
    y_zero    = (x > 0) & (y == 0)   # not in y dataset

    # Floor: one decade below the smallest nonzero value across both axes
    all_pos = pd.concat([x[x > 0], y[y > 0]])
    floor = 10 ** (np.floor(np.log10(all_pos.min())) - 0.5)

    # --- Main scatter (both > 0) ---
    x_ok = x[both_pos]
    y_ok = y[both_pos]
    lab_ok = labels[both_pos]
    cat_ok = categories[both_pos]

    order = sorted(set(cat_ok), key=lambda c: ZORDER.get(c, -1))
    for cat in order:
        sel = cat_ok == cat
        ax.scatter(
            x_ok[sel], y_ok[sel],
            c=PALETTE.get(cat, '#cccccc'),
            s=16, alpha=0.7, linewidths=0,
            label=cat, zorder=ZORDER.get(cat, 1) + 2,
            rasterized=True,
        )

    # --- Off-axis: x=0 (left-pointing triangles at x floor) ---
    if x_zero.any():
        x_xz = x_ok.iloc[:0]   # empty — marker x position is floor
        y_xz = y[x_zero]
        cat_xz = categories[x_zero]
        lab_xz = labels[x_zero]
        marker_x = floor * 3
        for cat in sorted(set(cat_xz), key=lambda c: ZORDER.get(c, -1)):
            sel = cat_xz == cat
            ax.scatter(
                [marker_x] * sel.sum(), y_xz[sel],
                marker='<', c=PALETTE.get(cat, '#cccccc'),
                s=45, linewidths=0.6, edgecolors='#333333',
                zorder=ZORDER.get(cat, 1) + 2,
            )
        for i in range(len(lab_xz)):
            ax.annotate(
                f"{lab_xz.iloc[i]} (x=0)",
                xy=(marker_x, y_xz.iloc[i]),
                xytext=(-4, 0), textcoords='offset points',
                fontsize=6, ha='right', va='center', color='#444444',
                path_effects=[pe.withStroke(linewidth=1.5, foreground='white')],
            )

    # --- Off-axis: y=0 (downward triangles at y floor) ---
    if y_zero.any():
        x_yz = x[y_zero]
        cat_yz = categories[y_zero]
        lab_yz = labels[y_zero]
        marker_y = floor * 3
        for cat in sorted(set(cat_yz), key=lambda c: ZORDER.get(c, -1)):
            sel = cat_yz == cat
            ax.scatter(
                x_yz[sel], [marker_y] * sel.sum(),
                marker='v', c=PALETTE.get(cat, '#cccccc'),
                s=45, linewidths=0.6, edgecolors='#333333',
                zorder=ZORDER.get(cat, 1) + 2,
            )
        for i in range(len(lab_yz)):
            ax.annotate(
                f"{lab_yz.iloc[i]} (y=0)",
                xy=(x_yz.iloc[i], marker_y),
                xytext=(4, -8), textcoords='offset points',
                fontsize=6, va='top', color='#444444',
                path_effects=[pe.withStroke(linewidth=1.5, foreground='white')],
            )

    # --- Axes and scales ---
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(left=floor)
    ax.set_ylim(bottom=floor)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12)

    # Identity line (y = x)
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(x_ok.max(), y_ok.max()) * 1.5,
    ]
    ax.plot(lims, lims, color='#888888', lw=1.0, ls='-', zorder=1, label='_identity')

    # Pearson r
    lx, ly = np.log10(x_ok), np.log10(y_ok)
    r, _ = stats.pearsonr(lx, ly)
    ax.text(
        0.04, 0.97, f'r = {r:.3f}',
        transform=ax.transAxes, va='top', fontsize=10,
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.8),
    )

    # Annotate non-other genes by |log2 fold-change from identity|
    non_other = cat_ok != 'other'
    x_ann = x_ok[non_other]
    y_ann = y_ok[non_other]
    lab_ann = lab_ok[non_other]
    log2fc = np.log2(y_ann.values) - np.log2(x_ann.values)
    n_each = n_annotate // 2
    top_idx = np.argsort(log2fc)[-n_each:]
    bot_idx = np.argsort(log2fc)[:n_each]
    for i in np.concatenate([top_idx, bot_idx]):
        ax.annotate(
            lab_ann.iloc[i],
            xy=(x_ann.iloc[i], y_ann.iloc[i]),
            xytext=(4, 0), textcoords='offset points',
            fontsize=6.5, va='center', color='#222222',
            path_effects=[pe.withStroke(linewidth=1.5, foreground='white')],
        )

    return r


# ---------------------------------------------------------------------------
# Concordance report
# ---------------------------------------------------------------------------

# All annotation columns from gene_metadata (after gene_id / gene_symbol)
_ANNOTATION_COLS = [
    'is_essential', 'aa_pway_enzyme', 'is_tf',
    'is_ribosomal_translation_mach', 'is_rnap_transcription_mach',
    'has_expression_adjustment', 'expression_adjustment_factor',
    'has_deg_rate_adjustment', 'deg_rate_adjustment_factor',
]


def _compute_log2fc(tpm_x: pd.Series, tpm_y: pd.Series) -> np.ndarray:
    """
    Compute log2(y / x) element-wise.

    Special cases:
      x > 0, y > 0  ->  normal log2 ratio
      x > 0, y = 0  ->  -inf  (gene silenced in experiment)
      x = 0, y > 0  ->  +inf  (gene absent from reference)
      x = 0, y = 0  ->  NaN   (no data in either dataset; excluded upstream)
    """
    x = tpm_x.values.astype(float)
    y = tpm_y.values.astype(float)
    with np.errstate(divide='ignore', invalid='ignore'):
        fc = np.where(
            (x == 0) & (y == 0), np.nan,
            np.where(x == 0, np.inf,
            np.where(y == 0, -np.inf,
            np.log2(y / x))),
        )
    return fc


def correlation_summary(
    merged: pd.DataFrame,
    flag_thresholds: tuple = (-1.0, -2.0),
) -> pd.DataFrame:
    """
    Table 1: Pearson/Spearman r and downregulation counts for gene subsets.

    Correlations are computed on log10(TPM) for genes where tpm_x > 0 AND
    tpm_y > 0 (n_corr). Downregulation counts (n_log2fc_lt_*) cover all genes
    in the subset where tpm_x > 0, including genes absent from the experiment
    (tpm_y = 0, which have log2fc = -inf and are always flagged).

    Parameters
    ----------
    merged : DataFrame containing tpm_x, tpm_y, and all gene_metadata columns.
    flag_thresholds : log2fc cutoffs for the count columns, e.g. (-1.0, -2.0)
                      produces n_log2fc_lt_1 and n_log2fc_lt_2. Values must be
                      negative. Default: (-1.0, -2.0).

    Returns
    -------
    DataFrame with columns: subset, n_total, n_corr, pearson_r, spearman_r,
    n_log2fc_lt_<t> for each threshold t.
    """
    df = merged.copy()
    df['_log2fc'] = _compute_log2fc(df['tpm_x'], df['tpm_y'])

    def _bool(col):
        return df[col].fillna(False).astype(bool)

    subsets = {
        'all in-model':           df,
        'essential':              df[_bool('is_essential')],
        'expression adjusted':    df[_bool('has_expression_adjustment')],
        'ribosome / translation': df[_bool('is_ribosomal_translation_mach')],
        'TF':                     df[_bool('is_tf')],
        'AA pathway enzyme':      df[_bool('aa_pway_enzyme')],
        'RNAP / sigma':           df[_bool('is_rnap_transcription_mach')],
        'tpm_ref > 100':          df[df['tpm_x'] > 100],
    }

    # Column name for each threshold: -1.0 -> "n_log2fc_lt_1", -1.5 -> "n_log2fc_lt_1.5"
    def _thresh_col(t):
        label = f'{abs(t):.1f}'.rstrip('0').rstrip('.')
        return f'n_log2fc_lt_{label}'

    rows = []
    for name, sub_all in subsets.items():
        both_pos = (sub_all['tpm_x'] > 0) & (sub_all['tpm_y'] > 0)
        sub = sub_all[both_pos]
        n_total = len(sub_all)
        n_corr  = len(sub)

        if n_corr < 3:
            pearson_r = spearman_r = np.nan
        else:
            lx = np.log10(sub['tpm_x'])
            ly = np.log10(sub['tpm_y'])
            pearson_r, _  = stats.pearsonr(lx, ly)
            spearman_r, _ = stats.spearmanr(lx, ly)

        row = {
            'subset':     name,
            'n_total':    n_total,
            'n_corr':     n_corr,
            'pearson_r':  round(pearson_r, 4) if not np.isnan(pearson_r) else np.nan,
            'spearman_r': round(spearman_r, 4) if not np.isnan(spearman_r) else np.nan,
        }
        for t in flag_thresholds:
            # Count genes in reference (tpm_x > 0) that are below threshold in expt,
            # including absent genes (tpm_y = 0, log2fc = -inf)
            eligible = sub_all[sub_all['tpm_x'] > 0]
            row[_thresh_col(t)] = int((eligible['_log2fc'] < t).sum())

        rows.append(row)
    return pd.DataFrame(rows)


def gene_log2fc_table(merged: pd.DataFrame) -> pd.DataFrame:
    """
    Full log2FC table for all in-model genes (unfiltered).

    Returns log2(tpm_y / tpm_x) for every gene, including:
      - Genes absent from experiment (tpm_y = 0): log2fc = -inf
      - Genes absent from reference (tpm_x = 0): log2fc = +inf
      - Genes absent from both: excluded upstream

    Sorted by category priority (expression-adjusted first, then essential,
    then other), then log2fc ascending within each group.

    Use this table for distributional analysis (e.g. violin plots). For the
    actionable flagged subset, use problematic_genes() instead.

    Parameters
    ----------
    merged : DataFrame with tpm_x, tpm_y, tpm_std_x, tpm_std_y, category, and
             all gene_metadata columns.

    Returns
    -------
    DataFrame with columns: gene_id, gene_symbol, category, log2fc,
    tpm_mean_ref, tpm_std_ref, tpm_mean_expt, tpm_std_expt,
    then all annotation columns from gene_metadata.
    """
    df = merged.copy()
    df['log2fc'] = _compute_log2fc(df['tpm_x'], df['tpm_y'])

    adj = df['has_expression_adjustment'].fillna(False).astype(bool)
    ess = df['is_essential'].fillna(False).astype(bool)
    df['_cat_priority'] = np.where(adj, 0, np.where(ess, 1, 2))
    df = (
        df
        .sort_values(['_cat_priority', 'log2fc'], ascending=[True, True])
        .drop(columns=['_cat_priority'])
    )

    df = df.rename(columns={
        'tpm_x':     'tpm_mean_ref',
        'tpm_std_x': 'tpm_std_ref',
        'tpm_y':     'tpm_mean_expt',
        'tpm_std_y': 'tpm_std_expt',
    })

    out_cols = (
        ['gene_id', 'gene_symbol', 'category', 'log2fc',
         'tpm_mean_ref', 'tpm_std_ref', 'tpm_mean_expt', 'tpm_std_expt']
        + [c for c in _ANNOTATION_COLS if c in df.columns]
    )
    return df[out_cols].reset_index(drop=True)


def problematic_genes(
    merged: pd.DataFrame,
    log2fc_threshold: float = -1.0,
) -> pd.DataFrame:
    """
    Table 2: In-model genes with lower expression in the experiment than in the
    reference (one-sided: we are looking for potential failure-causing dropouts).

    Genes are flagged when log2(tpm_y / tpm_x) < log2fc_threshold. Genes with
    tpm_y = 0 have log2fc = -inf and are always included. Genes absent from
    both datasets are excluded upstream.

    Rows are sorted by:
      1. Category priority: essential+adjusted > expression adjusted > essential > other
      2. log2fc ascending (most downregulated first; -inf rows lead each category)

    Parameters
    ----------
    merged : DataFrame with tpm_x, tpm_y, tpm_std_x, tpm_std_y, category, and
             all gene_metadata columns.
    log2fc_threshold : Flag genes with log2fc strictly below this value.
                       Default -1.0 (≥2-fold lower in experiment).

    Returns
    -------
    DataFrame with columns: gene_id, gene_symbol, category, log2fc,
    tpm_mean_ref, tpm_std_ref, tpm_mean_expt, tpm_std_expt,
    then all annotation columns from gene_metadata.
    """
    full = gene_log2fc_table(merged)
    return full[full['log2fc'] < log2fc_threshold].reset_index(drop=True)


def generate_concordance_report(
    x_id: str,
    y_id: str,
    log2fc_threshold: float = -1.0,
    flag_thresholds: tuple = (-1.0, -2.0),
    out_dir: Path = None,
) -> tuple:
    """
    Generate a concordance report comparing experiment dataset y_id to reference x_id.

    Writes two TSV files to out_dir (default: analysis/results/):
      <date>_concordance_corr__{x_id}__{y_id}.tsv   — Table 1 (correlations by subset)
      <date>_concordance_genes__{x_id}__{y_id}.tsv  — Table 2 (flagged genes)

    Also prints a summary to stdout.

    Parameters
    ----------
    x_id : Dataset ID for the reference (x-axis).
           Default CLI value: 'vecoli_m9_glucose_minus_aas'.
    y_id : Dataset ID for the experiment (y-axis).
    log2fc_threshold : Genes with log2fc < this value are flagged in Table 2.
                       Default -1.0 (≥2-fold lower in experiment).
    flag_thresholds : log2fc cutoffs used to compute downregulation count columns
                      in Table 1. Default (-1.0, -2.0).
    out_dir : Directory to write TSVs. Defaults to analysis/results/.

    Returns
    -------
    (table1, table2) as DataFrames.
    """
    if out_dir is None:
        out_dir = REPORTS_DIR

    x_df = pd.read_csv(DATA_DIR / f'{x_id}.tsv', sep='\t')
    y_df = pd.read_csv(DATA_DIR / f'{y_id}.tsv', sep='\t')
    meta = pd.read_csv(RESULTS_DIR / 'gene_metadata.tsv', sep='\t')

    merged = (
        x_df.rename(columns={'tpm_mean': 'tpm_x', 'tpm_std': 'tpm_std_x'})
        .merge(
            y_df.rename(columns={'tpm_mean': 'tpm_y', 'tpm_std': 'tpm_std_y'}),
            on='gene_id', how='inner',
        )
        .merge(meta, on='gene_id', how='left')
    )
    merged['tpm_x'] = merged['tpm_x'].fillna(0)
    merged['tpm_y'] = merged['tpm_y'].fillna(0)
    merged['category'] = merged.apply(assign_category, axis=1)

    in_model = merged['in_model'].fillna(False)
    merged = merged[in_model & ((merged['tpm_x'] > 0) | (merged['tpm_y'] > 0))].reset_index(drop=True)

    t1 = correlation_summary(merged, flag_thresholds=flag_thresholds)
    t2 = problematic_genes(merged, log2fc_threshold=log2fc_threshold)
    t3 = gene_log2fc_table(merged)

    date_prefix = date.today().strftime('%d%b%Y').lower()
    t1_path = out_dir / f'{date_prefix}_concordance_corr__{x_id}__{y_id}.tsv'
    t2_path = out_dir / f'{date_prefix}_concordance_genes__{x_id}__{y_id}.tsv'
    t3_path = out_dir / f'{date_prefix}_concordance_allgenes__{x_id}__{y_id}.tsv'
    t1.to_csv(t1_path, sep='\t', index=False)
    t2.to_csv(t2_path, sep='\t', index=False)
    t3.to_csv(t3_path, sep='\t', index=False)

    print(f"\n=== Concordance report: {y_id} vs {x_id} ===")
    print(f"\nTable 1 — Correlation by gene subset (log10 scale, both > 0):")
    print(t1.to_string(index=False))
    print(f"\nTable 2 — Flagged genes (log2fc < {log2fc_threshold}): {len(t2)} total")
    n_adj = t2['has_expression_adjustment'].fillna(False).astype(bool).sum()
    n_ess = (t2['is_essential'].fillna(False).astype(bool) & ~t2['has_expression_adjustment'].fillna(False).astype(bool)).sum()
    n_other = len(t2) - n_adj - n_ess
    for label, n in [('expression adjusted (any)', n_adj), ('essential (non-adjusted)', n_ess), ('other', n_other)]:
        if n:
            print(f"  {label}: {n}")
    print(f"\nSaved: {t1_path.name}")
    print(f"Saved: {t2_path.name}")
    print(f"Saved: {t3_path.name}")

    return t1, t2


# Compact column-name slugs for flattening Table 1 to a wide row
_SUBSET_SLUGS = {
    'all in-model':           'all',
    'essential':              'essential',
    'expression adjusted':    'expr_adj',
    'ribosome / translation': 'ribosome',
    'TF':                     'tf',
    'AA pathway enzyme':      'aa_enzyme',
    'RNAP / sigma':           'rnap',
    'tpm_ref > 100':          'tpm_gt100',
}


def flatten_summary(t1: pd.DataFrame, dataset_id: str) -> pd.DataFrame:
    """
    Pivot the 8-row correlation summary (Table 1) into a single wide row.

    Column names follow the pattern {subset_slug}__{metric}, e.g.
    essential__pearson_r, expr_adj__n_log2fc_lt_1. This format is suitable
    for stacking across multiple datasets to build a meta-analysis table.

    Parameters
    ----------
    t1 : DataFrame output of correlation_summary().
    dataset_id : Dataset identifier prepended as the first column.

    Returns
    -------
    Single-row DataFrame with dataset_id plus one column per (subset × metric).
    """
    metric_cols = [c for c in t1.columns if c != 'subset']
    row = {'dataset_id': dataset_id}
    for _, r in t1.iterrows():
        slug = _SUBSET_SLUGS.get(
            r['subset'],
            r['subset'].replace(' ', '_').replace('/', '_'),
        )
        for col in metric_cols:
            row[f'{slug}__{col}'] = r[col]
    return pd.DataFrame([row])


def generate_meta_summary(
    dataset_ids: list,
    x_id: str = 'vecoli_m9_glucose_minus_aas',
    parca_csvs: list = None,
    flag_thresholds: tuple = (-1.0, -2.0),
    out_dir: Path = None,
) -> pd.DataFrame:
    """
    Generate concordance reports for a list of datasets and combine them into
    a single wide meta-summary table for comparison against modeling outcomes.

    For each dataset_id, writes individual concordance TSVs to out_dir and
    accumulates a flattened Table 1 row. If parca_csvs are provided, joins
    parca_status, parca_error, and n_genes_filled_from_ref by dataset_id,
    placing outcome columns immediately after dataset_id.

    Parameters
    ----------
    dataset_ids : Dataset IDs to process. The reference (x_id) is skipped if
                  included.
    x_id : Reference dataset ID used as the x-axis in all comparisons.
           Default: 'vecoli_m9_glucose_minus_aas'.
    parca_csvs : Paths to parca summary CSVs whose modeling outcomes should be
                 joined. Multiple CSVs are concatenated; if a dataset_id appears
                 in more than one file, the last occurrence is kept.
    flag_thresholds : log2fc cutoffs for downregulation count columns.
                      Default: (-1.0, -2.0).
    out_dir : Output directory for individual concordance TSVs and the
              meta-summary TSV. Defaults to analysis/results/.

    Returns
    -------
    Wide DataFrame with one row per dataset.
    """
    if out_dir is None:
        out_dir = REPORTS_DIR

    rows = []
    for y_id in dataset_ids:
        if y_id == x_id:
            print(f"Skipping reference dataset: {y_id}")
            continue
        print(f"\nProcessing: {y_id}")
        t1, _ = generate_concordance_report(
            x_id=x_id,
            y_id=y_id,
            flag_thresholds=flag_thresholds,
            out_dir=out_dir,
        )
        rows.append(flatten_summary(t1, y_id))

    meta = pd.concat(rows, ignore_index=True)

    if parca_csvs:
        parca = pd.concat(
            [pd.read_csv(p) for p in parca_csvs], ignore_index=True
        ).drop_duplicates('dataset_id', keep='last')
        join_cols = [c for c in
                     ['dataset_id', 'parca_status', 'parca_error',
                      'sim_errors', 'n_genes_filled_from_ref']
                     if c in parca.columns]
        meta = meta.merge(parca[join_cols], on='dataset_id', how='left')

        # Derive a 3-level outcome:
        #   success      — parca completed AND no sim errors
        #   sim_failed   — parca completed BUT sim errors present
        #   parca_failed — parca did not complete
        parca_ok = meta['parca_status'] == 'COMPLETED'
        sim_ok   = meta['sim_errors'].isna() | (meta['sim_errors'].astype(str).str.strip() == '')
        meta['outcome'] = 'parca_failed'
        meta.loc[parca_ok & ~sim_ok,  'outcome'] = 'sim_failed'
        meta.loc[parca_ok &  sim_ok,  'outcome'] = 'success'

        # Place outcome columns immediately after dataset_id
        outcome_cols = ['outcome'] + [c for c in join_cols if c != 'dataset_id']
        feature_cols = [c for c in meta.columns if c not in ['dataset_id'] + outcome_cols]
        meta = meta[['dataset_id'] + outcome_cols + feature_cols]

    date_prefix = date.today().strftime('%d%b%Y').lower()
    out_path = out_dir / f'{date_prefix}_meta_summary.tsv'
    meta.to_csv(out_path, sep='\t', index=False)
    print(f"\nMeta-summary ({len(meta)} datasets) saved: {out_path.name}")

    return meta


# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Compare TPM of two datasets.')
    parser.add_argument('--x', default='vecoli_m9_glucose_minus_aas',
                        help='Dataset ID for x-axis / reference (default: vecoli reference)')
    parser.add_argument('--y', default='gbw_vegas_wt_m9glc_34h_no_ssrA',
                        help='Dataset ID for y-axis / experiment')
    parser.add_argument('--n_annotate', type=int, default=16,
                        help='Number of outlier genes to label (split evenly above/below identity)')
    parser.add_argument('--report', action='store_true',
                        help='Generate concordance report TSVs in addition to the scatter plot')
    parser.add_argument('--log2fc_threshold', type=float, default=-1.0,
                        help='log2fc cutoff for flagging downregulated genes in the report '
                             '(default: -1.0, i.e. ≥2-fold lower in experiment)')
    args = parser.parse_args()

    # Load data
    x_df = pd.read_csv(DATA_DIR / f'{args.x}.tsv', sep='\t')
    y_df = pd.read_csv(DATA_DIR / f'{args.y}.tsv', sep='\t')
    meta = pd.read_csv(RESULTS_DIR / 'gene_metadata.tsv', sep='\t')

    merged = (
        x_df.rename(columns={'tpm_mean': 'tpm_x', 'tpm_std': 'tpm_std_x'})
        .merge(
            y_df.rename(columns={'tpm_mean': 'tpm_y', 'tpm_std': 'tpm_std_y'}),
            on='gene_id', how='inner',
        )
        .merge(meta[['gene_id', 'gene_symbol', 'in_model', 'is_essential', 'has_expression_adjustment']],
               on='gene_id', how='left')
    )
    merged['tpm_x'] = merged['tpm_x'].fillna(0)
    merged['tpm_y'] = merged['tpm_y'].fillna(0)
    merged['category'] = merged.apply(assign_category, axis=1)

    # Restrict to in-model genes; drop both-zero
    in_model = merged['in_model'].fillna(False)
    merged = merged[in_model & ((merged['tpm_x'] > 0) | (merged['tpm_y'] > 0))].reset_index(drop=True)

    print(f"x dataset ({args.x}): {len(x_df)} genes")
    print(f"y dataset ({args.y}): {len(y_df)} genes")
    print(f"In-model genes with data in at least one dataset: {len(merged)}")
    print(f"  both > 0:  {((merged['tpm_x']>0) & (merged['tpm_y']>0)).sum()}")
    print(f"  x=0, y>0:  {((merged['tpm_x']==0) & (merged['tpm_y']>0)).sum()}")
    print(f"  x>0, y=0:  {((merged['tpm_x']>0) & (merged['tpm_y']==0)).sum()}")

    fig, ax = plt.subplots(figsize=(7, 7))

    r = tpm_comparison_scatter(
        ax,
        x=merged['tpm_x'],
        y=merged['tpm_y'],
        labels=merged['gene_symbol'],
        categories=merged['category'],
        xlabel=f'TPM — {args.x}',
        ylabel=f'TPM — {args.y}',
        title=f'Dataset TPM comparison\n{args.x}  vs  {args.y}',
        n_annotate=args.n_annotate,
    )

    # Legend
    present = merged['category'].unique()
    handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=PALETTE[c],
                   markersize=7, label=c)
        for c in PALETTE if c in present
    ]
    handles += [
        plt.Line2D([0], [0], color='#888888', lw=1.0, label='identity (y=x)'),
    ]
    ax.legend(handles=handles, fontsize=8, loc='upper left', framealpha=0.9)

    fig.tight_layout()

    out_name = f'tpm_compare__{args.x}__{args.y}.png'
    fig.savefig(FIGURES_DIR / out_name, dpi=150)
    plt.close(fig)
    print(f"\nSaved: figures/{out_name}")
    print(f"Pearson r (log-log, both>0): {r:.3f}")

    if args.report:
        generate_concordance_report(
            x_id=args.x,
            y_id=args.y,
            log2fc_threshold=args.log2fc_threshold,
        )


if __name__ == '__main__':
    main()
