"""
5_tpm_to_model_expression.py

Scatter plots of TPM (vecoli_m9_glucose_minus_aas reference dataset) vs
model-fitted expression quantities from parca:
  - TPM vs mRNA steady-state mean count  (mrna_ss_mean_count)
  - TPM vs monomer steady-state mean count  (monomer_ss_mean_count)

Run from analysis/ directory:
    python 5_tpm_to_model_expression.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from pathlib import Path
from scipy import stats

HERE = Path(__file__).parent
DATA_DIR = HERE / '../data_formatted'
RESULTS_DIR = HERE / 'model_results'
FIGURES_DIR = HERE / 'figures'
FIGURES_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Load and merge
# ---------------------------------------------------------------------------

tpm = pd.read_csv(DATA_DIR / 'vecoli_m9_glucose_minus_aas.tsv', sep='\t')
abund = pd.read_csv(RESULTS_DIR / 'ref_basal_abundances.tsv', sep='\t')
meta = pd.read_csv(RESULTS_DIR / 'gene_metadata.tsv', sep='\t')

df = (
    tpm
    .merge(abund[['gene_id', 'gene_symbol', 'mrna_ss_mean_count', 'monomer_ss_mean_count']], on='gene_id', how='inner')
    .merge(meta.drop(columns='gene_symbol'), on='gene_id', how='left')
)

# Only genes that are in the model and have finite values on both axes
in_model = df['in_model'].fillna(False)
df_plot = df[in_model].copy()

print(f"Total genes in ref dataset:  {len(tpm)}")
print(f"Genes in model:              {in_model.sum()}")
print(f"Genes with mRNA count:       {df_plot['mrna_ss_mean_count'].notna().sum()}")
print(f"Genes with monomer count:    {df_plot['monomer_ss_mean_count'].notna().sum()}")

# ---------------------------------------------------------------------------
# Gene-category coloring
# Priority order: adjustment > rnap > ribosome > tf > aa_enzyme > plain
# ---------------------------------------------------------------------------

def assign_category(row):
    if row.get('has_expression_adjustment', False):
        return 'expression adjusted'
    if row.get('is_rnap_transcription_mach', False):
        return 'RNAP / sigma'
    if row.get('is_ribosomal_translation_mach', False):
        return 'ribosome / translation'
    if row.get('is_tf', False):
        return 'TF'
    if row.get('aa_pway_enzyme', False):
        return 'AA pathway enzyme'
    return 'other'

df_plot['category'] = df_plot.apply(assign_category, axis=1)

PALETTE = {
    'other':                   '#aaaaaa',
    'TF':                      '#4c78a8',
    'AA pathway enzyme':       '#f58518',
    'ribosome / translation':  '#54a24b',
    'RNAP / sigma':            '#b279a2',
    'expression adjusted':     '#e45756',
}
ZORDER = {cat: i for i, cat in enumerate(PALETTE)}  # draw special cats on top

# ---------------------------------------------------------------------------
# Helper: one scatter panel
# ---------------------------------------------------------------------------

def scatter_panel(
    ax, x, y, labels, categories,
    xlabel, ylabel, title,
    n_annotate=15,
):
    """
    Log-log scatter colored by category.
    Points more than 2 decades below the 1st-percentile floor are clipped:
    drawn as downward triangles at the floor with their actual value annotated.
    Annotates the n_annotate most extreme residual genes (above + below fit line)
    among non-clipped points.
    """
    mask = np.isfinite(np.log10(x)) & np.isfinite(np.log10(y)) & (x > 0) & (y > 0)
    x = x[mask].reset_index(drop=True)
    y = y[mask].reset_index(drop=True)
    labels = labels[mask].reset_index(drop=True)
    categories = categories[mask].reset_index(drop=True)

    # Determine clip floor: one decade below the 1st-percentile value
    floor = 10 ** (np.floor(np.log10(np.percentile(y, 1))) - 1)
    clipped = y < floor

    x_ok, y_ok = x[~clipped], y[~clipped]
    lab_ok, cat_ok = labels[~clipped], categories[~clipped]
    x_cl, y_cl = x[clipped], y[clipped]
    lab_cl, cat_cl = labels[clipped], categories[clipped]

    # --- Draw in-range points ---
    order = sorted(set(cat_ok), key=lambda c: ZORDER.get(c, -1))
    for cat in order:
        sel = cat_ok == cat
        ax.scatter(
            x_ok[sel], y_ok[sel],
            c=PALETTE.get(cat, '#aaaaaa'),
            s=18, alpha=0.7, linewidths=0,
            label=cat, zorder=ZORDER.get(cat, 1) + 2,
            rasterized=True,
        )

    # --- Draw clipped points as downward triangles just above the floor ---
    marker_y = floor * 3   # sits visually just above the bottom axis
    for cat in sorted(set(cat_cl), key=lambda c: ZORDER.get(c, -1)):
        sel = cat_cl == cat
        ax.scatter(
            x_cl[sel], [marker_y] * sel.sum(),
            marker='v', c=PALETTE.get(cat, '#aaaaaa'),
            s=50, linewidths=0.6, edgecolors='#333333',
            zorder=ZORDER.get(cat, 1) + 2,
        )
    for i in range(len(x_cl)):
        ax.annotate(
            f"{lab_cl.iloc[i]} ({y_cl.iloc[i]:.1e})",
            xy=(x_cl.iloc[i], marker_y),
            xytext=(4, -10), textcoords='offset points',
            fontsize=6.5, va='top', color='#222222',
            path_effects=[pe.withStroke(linewidth=1.5, foreground='white')],
        )

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(bottom=floor)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12)

    # Pearson r on log-transformed in-range values
    lx, ly = np.log10(x_ok), np.log10(y_ok)
    r, _ = stats.pearsonr(lx, ly)
    ax.text(
        0.04, 0.97, f'r = {r:.3f}',
        transform=ax.transAxes, va='top', fontsize=10,
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.8),
    )

    # OLS fit line (on in-range points)
    slope, intercept, *_ = stats.linregress(lx, ly)
    xlim = np.array([x_ok.min() * 0.8, x_ok.max() * 1.2])
    ax.plot(
        xlim,
        10 ** (slope * np.log10(xlim) + intercept),
        color='#333333', lw=1.2, ls='--', zorder=1,
    )

    # Annotate residual outliers among in-range points
    residuals = ly - (slope * lx + intercept)
    n_each = n_annotate // 2
    top_idx = np.argsort(residuals.values)[-n_each:]
    bot_idx = np.argsort(residuals.values)[:n_each]
    for i in np.concatenate([top_idx, bot_idx]):
        ax.annotate(
            lab_ok.iloc[i],
            xy=(x_ok.iloc[i], y_ok.iloc[i]),
            xytext=(4, 0), textcoords='offset points',
            fontsize=6.5, va='center', color='#222222',
            path_effects=[pe.withStroke(linewidth=1.5, foreground='white')],
        )

    return r


# ---------------------------------------------------------------------------
# Figure 1: TPM vs mRNA ss mean count
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(7, 6))

mrna_mask = df_plot['mrna_ss_mean_count'].notna() & (df_plot['mrna_ss_mean_count'] > 0) & (df_plot['tpm_mean'] > 0)
sub = df_plot[mrna_mask].reset_index(drop=True)

r_mrna = scatter_panel(
    ax,
    x=sub['tpm_mean'],
    y=sub['mrna_ss_mean_count'],
    labels=sub['gene_symbol'],
    categories=sub['category'],
    xlabel='TPM (ref dataset, mean across replicates)',
    ylabel='mRNA steady-state mean count (parca fitted)',
    title='TPM → mRNA count\n(vecoli_m9_glucose_minus_aas, in-model genes)',
)

# Legend — only categories actually present
present = sub['category'].unique()
handles = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=PALETTE[c],
               markersize=7, label=c)
    for c in PALETTE if c in present
]
ax.legend(handles=handles, fontsize=8, loc='lower right', framealpha=0.9)
fig.tight_layout()
fig.savefig(FIGURES_DIR / 'tpm_vs_mrna_count.png', dpi=150)
plt.close(fig)
print("Saved: figures/tpm_vs_mrna_count.png")

# ---------------------------------------------------------------------------
# Figure 2: TPM vs monomer ss mean count
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(7, 6))

mono_mask = df_plot['monomer_ss_mean_count'].notna() & (df_plot['monomer_ss_mean_count'] > 0) & (df_plot['tpm_mean'] > 0)
sub2 = df_plot[mono_mask].reset_index(drop=True)

r_mono = scatter_panel(
    ax,
    x=sub2['tpm_mean'],
    y=sub2['monomer_ss_mean_count'],
    labels=sub2['gene_symbol'],
    categories=sub2['category'],
    xlabel='TPM (ref dataset, mean across replicates)',
    ylabel='Monomer steady-state mean count (parca fitted)',
    title='TPM → Protein count\n(vecoli_m9_glucose_minus_aas, in-model genes)',
)

present2 = sub2['category'].unique()
handles2 = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=PALETTE[c],
               markersize=7, label=c)
    for c in PALETTE if c in present2
]
ax.legend(handles=handles2, fontsize=8, loc='lower right', framealpha=0.9)
fig.tight_layout()
fig.savefig(FIGURES_DIR / 'tpm_vs_monomer_count.png', dpi=150)
plt.close(fig)
print("Saved: figures/tpm_vs_monomer_count.png")

# ---------------------------------------------------------------------------
# Figure 3: side-by-side summary (convenient single-image version)
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, (ymask_col, ylabel, title_suffix) in zip(axes, [
    ('mrna_ss_mean_count',    'mRNA steady-state mean count',    'mRNA'),
    ('monomer_ss_mean_count', 'Monomer steady-state mean count', 'Protein'),
]):
    smask = df_plot[ymask_col].notna() & (df_plot[ymask_col] > 0) & (df_plot['tpm_mean'] > 0)
    sub_ax = df_plot[smask].reset_index(drop=True)
    scatter_panel(
        ax,
        x=sub_ax['tpm_mean'],
        y=sub_ax[ymask_col],
        labels=sub_ax['gene_symbol'],
        categories=sub_ax['category'],
        xlabel='TPM (ref dataset)',
        ylabel=ylabel,
        title=f'TPM → {title_suffix} count',
        n_annotate=10,
    )

# Shared legend on the right panel
present_all = df_plot['category'].unique()
handles_all = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=PALETTE[c],
               markersize=7, label=c)
    for c in PALETTE if c in present_all
]
axes[1].legend(handles=handles_all, fontsize=8, loc='lower right', framealpha=0.9)

fig.suptitle(
    'Reference dataset: TPM → model expression mapping\n'
    '(vecoli_m9_glucose_minus_aas, in-model genes only)',
    fontsize=12, y=1.01,
)
fig.tight_layout()
fig.savefig(FIGURES_DIR / 'tpm_vs_model_expression.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print("Saved: figures/tpm_vs_model_expression.png")

print(f"\nPearson r (log-log)  TPM vs mRNA count:    {r_mrna:.3f}")
print(f"Pearson r (log-log)  TPM vs protein count: {r_mono:.3f}")

# ---------------------------------------------------------------------------
# Figures 4–6: essential × expression-adjusted view
# ---------------------------------------------------------------------------

def assign_essential_category(row):
    ess = bool(row.get('is_essential', False))
    adj = bool(row.get('has_expression_adjustment', False))
    if ess and adj:
        return 'essential + adjusted'
    if ess:
        return 'essential'
    if adj:
        return 'expression adjusted'
    return 'other'

df_plot['ess_category'] = df_plot.apply(assign_essential_category, axis=1)

ESS_PALETTE = {
    'other':                  '#cccccc',
    'essential':              '#4c78a8',
    'expression adjusted':    '#e45756',
    'essential + adjusted':   '#f58518',
}
ESS_ZORDER = {cat: i for i, cat in enumerate(ESS_PALETTE)}

def scatter_panel_ess(ax, x, y, labels, categories, xlabel, ylabel, title, n_annotate=15):
    """Same as scatter_panel but uses ESS_PALETTE / ESS_ZORDER."""
    mask = np.isfinite(np.log10(x)) & np.isfinite(np.log10(y)) & (x > 0) & (y > 0)
    x = x[mask].reset_index(drop=True)
    y = y[mask].reset_index(drop=True)
    labels = labels[mask].reset_index(drop=True)
    categories = categories[mask].reset_index(drop=True)

    floor = 10 ** (np.floor(np.log10(np.percentile(y, 1))) - 1)
    clipped = y < floor

    x_ok, y_ok = x[~clipped], y[~clipped]
    lab_ok, cat_ok = labels[~clipped], categories[~clipped]
    x_cl, y_cl = x[clipped], y[clipped]
    lab_cl, cat_cl = labels[clipped], categories[clipped]

    order = sorted(set(cat_ok), key=lambda c: ESS_ZORDER.get(c, -1))
    for cat in order:
        sel = cat_ok == cat
        ax.scatter(
            x_ok[sel], y_ok[sel],
            c=ESS_PALETTE.get(cat, '#cccccc'),
            s=18, alpha=0.7, linewidths=0,
            label=cat, zorder=ESS_ZORDER.get(cat, 1) + 2,
            rasterized=True,
        )

    marker_y = floor * 3
    for cat in sorted(set(cat_cl), key=lambda c: ESS_ZORDER.get(c, -1)):
        sel = cat_cl == cat
        ax.scatter(
            x_cl[sel], [marker_y] * sel.sum(),
            marker='v', c=ESS_PALETTE.get(cat, '#cccccc'),
            s=50, linewidths=0.6, edgecolors='#333333',
            zorder=ESS_ZORDER.get(cat, 1) + 2,
        )
    for i in range(len(x_cl)):
        ax.annotate(
            f"{lab_cl.iloc[i]} ({y_cl.iloc[i]:.1e})",
            xy=(x_cl.iloc[i], marker_y),
            xytext=(4, -10), textcoords='offset points',
            fontsize=6.5, va='top', color='#222222',
            path_effects=[pe.withStroke(linewidth=1.5, foreground='white')],
        )

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(bottom=floor)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12)

    lx, ly = np.log10(x_ok), np.log10(y_ok)
    r, _ = stats.pearsonr(lx, ly)
    ax.text(
        0.04, 0.97, f'r = {r:.3f}',
        transform=ax.transAxes, va='top', fontsize=10,
        bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.8),
    )

    slope, intercept, *_ = stats.linregress(lx, ly)
    xlim = np.array([x_ok.min() * 0.8, x_ok.max() * 1.2])
    ax.plot(
        xlim,
        10 ** (slope * np.log10(xlim) + intercept),
        color='#333333', lw=1.2, ls='--', zorder=1,
    )

    residuals = ly - (slope * lx + intercept)
    n_each = n_annotate // 2
    top_idx = np.argsort(residuals.values)[-n_each:]
    bot_idx = np.argsort(residuals.values)[:n_each]
    for i in np.concatenate([top_idx, bot_idx]):
        ax.annotate(
            lab_ok.iloc[i],
            xy=(x_ok.iloc[i], y_ok.iloc[i]),
            xytext=(4, 0), textcoords='offset points',
            fontsize=6.5, va='center', color='#222222',
            path_effects=[pe.withStroke(linewidth=1.5, foreground='white')],
        )

    return r


def make_legend(ax, palette, present):
    handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=palette[c],
                   markersize=7, label=c)
        for c in palette if c in present
    ]
    ax.legend(handles=handles, fontsize=8, loc='lower right', framealpha=0.9)


# Figure 4: TPM vs mRNA — essential view
fig, ax = plt.subplots(figsize=(7, 6))
mrna_mask = df_plot['mrna_ss_mean_count'].notna() & (df_plot['mrna_ss_mean_count'] > 0) & (df_plot['tpm_mean'] > 0)
sub_e = df_plot[mrna_mask].reset_index(drop=True)
scatter_panel_ess(
    ax, sub_e['tpm_mean'], sub_e['mrna_ss_mean_count'],
    sub_e['gene_symbol'], sub_e['ess_category'],
    xlabel='TPM (ref dataset, mean across replicates)',
    ylabel='mRNA steady-state mean count (parca fitted)',
    title='TPM → mRNA count  [essential / adjusted]\n(vecoli_m9_glucose_minus_aas)',
)
make_legend(ax, ESS_PALETTE, sub_e['ess_category'].unique())
fig.tight_layout()
fig.savefig(FIGURES_DIR / 'tpm_vs_mrna_count_essential.png', dpi=150)
plt.close(fig)
print("Saved: figures/tpm_vs_mrna_count_essential.png")

# Figure 5: TPM vs monomer — essential view
fig, ax = plt.subplots(figsize=(7, 6))
mono_mask = df_plot['monomer_ss_mean_count'].notna() & (df_plot['monomer_ss_mean_count'] > 0) & (df_plot['tpm_mean'] > 0)
sub_e2 = df_plot[mono_mask].reset_index(drop=True)
scatter_panel_ess(
    ax, sub_e2['tpm_mean'], sub_e2['monomer_ss_mean_count'],
    sub_e2['gene_symbol'], sub_e2['ess_category'],
    xlabel='TPM (ref dataset, mean across replicates)',
    ylabel='Monomer steady-state mean count (parca fitted)',
    title='TPM → Protein count  [essential / adjusted]\n(vecoli_m9_glucose_minus_aas)',
)
make_legend(ax, ESS_PALETTE, sub_e2['ess_category'].unique())
fig.tight_layout()
fig.savefig(FIGURES_DIR / 'tpm_vs_monomer_count_essential.png', dpi=150)
plt.close(fig)
print("Saved: figures/tpm_vs_monomer_count_essential.png")

# Figure 6: side-by-side — essential view
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax, (ymask_col, ylabel, title_suffix) in zip(axes, [
    ('mrna_ss_mean_count',    'mRNA steady-state mean count',    'mRNA'),
    ('monomer_ss_mean_count', 'Monomer steady-state mean count', 'Protein'),
]):
    smask = df_plot[ymask_col].notna() & (df_plot[ymask_col] > 0) & (df_plot['tpm_mean'] > 0)
    sub_ax = df_plot[smask].reset_index(drop=True)
    scatter_panel_ess(
        ax, sub_ax['tpm_mean'], sub_ax[ymask_col],
        sub_ax['gene_symbol'], sub_ax['ess_category'],
        xlabel='TPM (ref dataset)',
        ylabel=ylabel,
        title=f'TPM → {title_suffix} count  [essential / adjusted]',
        n_annotate=10,
    )
make_legend(axes[1], ESS_PALETTE, df_plot['ess_category'].unique())
fig.suptitle(
    'Reference dataset: TPM → model expression  [essential × adjusted view]\n'
    '(vecoli_m9_glucose_minus_aas, in-model genes only)',
    fontsize=12, y=1.01,
)
fig.tight_layout()
fig.savefig(FIGURES_DIR / 'tpm_vs_model_expression_essential.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print("Saved: figures/tpm_vs_model_expression_essential.png")
