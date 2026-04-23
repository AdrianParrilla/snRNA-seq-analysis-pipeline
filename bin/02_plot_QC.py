#!/opt/env/bin/python

import os
import warnings
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation as mad


def get_qc_thresholds(df_genes, groupby_key)-> pd.DataFrame:

    def valid_vals(x):
        """Helper: drop NaN and 0 before computing stats."""
        return x[(~x.isna()) & (x != 0)]

    qc_thresholds_upper = df_genes.groupby(groupby_key).agg({
        'log1p_total_counts': lambda x: np.median(valid_vals(x)) + 5 * mad(valid_vals(x)),
        'log1p_n_genes_by_counts': lambda x: np.median(valid_vals(x)) + 5 * mad(valid_vals(x)),
        'pct_counts_in_top_20_genes': lambda x: np.median(valid_vals(x)) + 5 * mad(valid_vals(x)),
        'pct_counts_mt': lambda x: min(np.median(valid_vals(x)) + 3 * mad(valid_vals(x)), 5)
    })

    qc_thresholds_lower = df_genes.groupby(groupby_key).agg({
        'log1p_total_counts': lambda x: np.median(valid_vals(x)) - 5 * mad(valid_vals(x)),
        'log1p_n_genes_by_counts': lambda x: max(np.median(valid_vals(x)) - 5 * mad(valid_vals(x)), np.log1p(100)),
        'pct_counts_in_top_20_genes': lambda x: max(np.median(valid_vals(x)) - 5 * mad(valid_vals(x)), 0)
    })
    
    return qc_thresholds_upper, qc_thresholds_lower


def draw_vertical_thresholds(ax, qc_upper_series, qc_lower_series=None,
                             y_ticks=None, y_labels=None,
                             color='red', lw=1.5, height_frac=0.25):
    """
    Draw vertical segments at x positions given by qc_upper_series (and qc_lower_series).
    y_ticks / y_labels are arrays of category positions & names (from the reference axis).
    """
    # If ticks/labels not provided, try to read them from the axis (fallback)
    if y_ticks is None or y_labels is None:
        ax.figure.canvas.draw()  # ensure positions exist
        y_ticks = ax.get_yticks()
        y_labels = [t.get_text() for t in ax.get_yticklabels()]

    if len(y_ticks) == 0:
        return

    # compute half-height from tick spacing
    spacing = np.min(np.abs(np.diff(y_ticks))) if len(y_ticks) > 1 else 1.0
    half_h = spacing * height_frac

    # reindex qc series to follow the y_labels order; this avoids key errors
    qc_up = qc_upper_series.reindex(y_labels)
    qc_low = qc_lower_series.reindex(y_labels) if qc_lower_series is not None else None

    for y, samp in zip(y_ticks, y_labels):
        if samp == '' or samp not in qc_up.index:
            continue
        x_up = qc_up.loc[samp]
        if pd.notna(x_up):
            ax.plot([x_up, x_up], [y - half_h, y + half_h], color=color, lw=lw)
        if qc_low is not None:
            x_low = qc_low.loc[samp]
            if pd.notna(x_low):
                ax.plot([x_low, x_low], [y - half_h, y + half_h], color=color, lw=lw)


def get_colors(df, color_by):

    '''
    Get the number of colors needed from the color palette
    '''
    color_palette = ["#fd756f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db",
                      "#fdcce5", "#bba387", "#c0d9af", "#b5b5b5", "#606060"] 

    n_colors = df[color_by].nunique()

    color_palette = color_palette[0:n_colors]

    return color_palette


def save_thresholds(qc_thresholds_upper, qc_thresholds_lower):
    
    thresholds = pd.merge(qc_thresholds_upper, qc_thresholds_lower, on='sample', suffixes=('_upper', '_lower'), how='inner')

    # reorder columns
    thresholds = thresholds[['log1p_n_genes_by_counts_lower', 'log1p_total_counts_upper', 'pct_counts_in_top_20_genes_lower', 'pct_counts_in_top_20_genes_upper', 'pct_counts_mt']]

    thresholds.to_csv('QC_thresholds.csv', index=True)


def save_QC_summary(df_genes, color_by):

    df_genes['scDblFinder_class'] = pd.Categorical(df_genes['scDblFinder_class']).rename_categories({1: 'singlet', 2: 'doublet'})

    counts = df_genes.groupby(['sample', 'outlier', 'scDblFinder_class', color_by],observed= True ).size().reset_index(name='total')

    counts.to_csv('QC_summary.csv', index=True)



def plot_qc(adata, sample_key, color_by, figsize: tuple, custom_order: list = None):

    warnings.simplefilter(action='ignore', category=FutureWarning)

    df_genes = adata.obs.copy()

    if custom_order:

        df_genes[sample_key] = pd.Categorical(df_genes[sample_key], categories=custom_order, ordered=True)

        df_genes = df_genes.sort_values(sample_key)

    if figsize is None:
        n_samples = df_genes[sample_key].nunique()
        
        # Define how much vertical space each sample gets (in inches)
        height_per_sample = 0.3 
        base_padding = 0 # Padding for x-labels, top legend, and margins
        
        # Calculate dynamic height (ensure a minimum height of 5 inches)
        calc_height = max(5.0, (n_samples * height_per_sample) + base_padding)
        
        figsize = (14, calc_height)
    
    qc_thresholds_upper, qc_thresholds_lower = get_qc_thresholds(df_genes, groupby_key=sample_key)

    save_thresholds(qc_thresholds_upper, qc_thresholds_lower)
    save_QC_summary(df_genes, color_by)

    color_palette = get_colors(df_genes, color_by)

    fig, axes = plt.subplots(1, 5, figsize=figsize, sharey=True, gridspec_kw={'wspace': 0.1})

    ax = axes[0]
    sns.violinplot(data=df_genes, x="log1p_total_counts", y=sample_key, hue=color_by, ax=ax, inner_kws={'box_width':3}, palette=color_palette)

    ax.set_xlabel('Nº of UMIs (log1p)', size=14, labelpad=10)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_ylabel('')

    if ax.get_legend():
        ax.get_legend().remove()

    ax.grid(False)

    # ensure the canvas is drawn so get_yticks() is accurate
    fig.canvas.draw()
    ref_y_ticks = ax.get_yticks()
    ref_y_labels = [t.get_text() for t in ax.get_yticklabels()]

    # Draw thresholds on first axis using reference ticks
    draw_vertical_thresholds(ax,
        qc_upper_series=qc_thresholds_upper['log1p_total_counts'],
        qc_lower_series=qc_thresholds_lower['log1p_total_counts'],
        y_ticks=ref_y_ticks, y_labels=ref_y_labels)


    ax = axes[1]
    sns.violinplot(data=df_genes, x="log1p_n_genes_by_counts", y=sample_key, hue=color_by, ax=ax, inner_kws={'box_width':3}, palette=color_palette)

    draw_vertical_thresholds(ax,
        qc_upper_series=qc_thresholds_upper['log1p_n_genes_by_counts'],
        qc_lower_series=qc_thresholds_lower['log1p_n_genes_by_counts'],
        y_ticks=ref_y_ticks, y_labels=ref_y_labels)

    ax.set_xlabel('Nº of genes (log1p)', size=14, labelpad=10)
    if ax.get_legend():
        ax.get_legend().remove()
    ax.grid(False)

    ax = axes[2]
    sns.violinplot(data=df_genes, x="pct_counts_in_top_20_genes", y=sample_key, hue=color_by, ax=ax, inner_kws={'box_width':3}, palette=color_palette)

    draw_vertical_thresholds(ax,
        qc_upper_series=qc_thresholds_upper['pct_counts_in_top_20_genes'],
        qc_lower_series=qc_thresholds_lower['pct_counts_in_top_20_genes'],
        y_ticks=ref_y_ticks, y_labels=ref_y_labels)

    ax.set_xlabel('% counts top 20 genes', size=14, labelpad=10)
    ax.set_xlim(-5, 115)
    if ax.get_legend():
        ax.get_legend().remove()
    ax.grid(False)


    ax = axes[3]
    sns.violinplot(data=df_genes, x="pct_counts_mt", y=sample_key, hue=color_by, ax=ax, inner_kws={'box_width':3}, palette=color_palette)

    draw_vertical_thresholds(ax,
        qc_upper_series=qc_thresholds_upper['pct_counts_mt'],
        qc_lower_series=None,
        y_ticks=ref_y_ticks, y_labels=ref_y_labels)

    ax.set_xlabel('%mt counts', size=14, labelpad=10)
    ax.get_legend().remove()
    ax.grid(False)

    ax = axes[4]
    if color_by == sample_key:
        batch_samples = df_genes[[sample_key]].drop_duplicates().reset_index(drop=True)

    else:
        batch_samples = df_genes[[sample_key, color_by]].drop_duplicates().reset_index(drop=True)
        
    n_cells = df_genes.groupby([sample_key, 'outlier']).size().reset_index(name='count')
    n_cells = n_cells.merge(batch_samples, on=sample_key)

    total_n_cells = n_cells.groupby(sample_key)['count'].sum().reset_index()
    total_n_cells = total_n_cells.merge(batch_samples, on=sample_key)

    n_cell_outliers = n_cells[n_cells['outlier'] == True]

    sns.barplot(data=total_n_cells, x='count', y=sample_key, hue=color_by, orient='h', edgecolor='black', linewidth=0.5, palette=color_palette)
    sns.barplot(data=n_cell_outliers, x='count', y=sample_key, color="#BFBFBF", orient='h', edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Nº of cells', size=14, labelpad=10)
    if ax.get_legend():
        ax.get_legend().remove()
    ax.grid(False)


    handles, labels = axes[0].get_legend_handles_labels()

    if len(labels) > 0:
        fig.legend(handles, labels, loc='lower center',
                bbox_to_anchor=(0.5, 0.9),
                ncol=len(labels), frameon=False)

    sns.set_style('ticks')
    plt.subplots_adjust(hspace=0.02)
    
    fig.savefig("QC_violin_plots.png", dpi=300, bbox_inches='tight')
    print("QC violin plots generated!")


def doublet_count(adata, sample_key, color_by):

    doublet_count = adata.obs[[sample_key, 'scDblFinder_class']].copy()

    doublet_count['scDblFinder_class'] = pd.Categorical(doublet_count['scDblFinder_class']).rename_categories({1: 'singlet', 2: 'doublet'})

    doublet_count = doublet_count.groupby([sample_key, 'scDblFinder_class'], sort=False).size().reset_index(name='counts')

    total_counts = doublet_count.groupby(sample_key)['counts'].sum().reset_index(name='total')

    doublet_count = doublet_count.merge(total_counts, on=sample_key)
    doublet_count = doublet_count.merge(adata.obs[[sample_key, color_by]].drop_duplicates(), on=sample_key)

    return doublet_count


def plot_doublets(adata, sample_key, color_by: str, figsize:tuple = None):

    doublet_counts = doublet_count(adata, sample_key, color_by)

    if figsize is None:
        n_samples = adata.obs[sample_key].nunique()
        
        # Define how much horizontal space each sample gets (in inches)
        width_per_sample = 0.35
        base_padding = 1 # Padding for x-labels, top legend, and margins
        
        # Calculate dynamic width (ensure a minimum width of 3 inches)
        calc_width = max(3.0, (n_samples * width_per_sample) + base_padding)
        
        figsize = (calc_width, 4)

    fig = plt.figure(figsize=figsize)

    color_palette = get_colors(doublet_counts, color_by)

    sns.barplot(data=doublet_counts, x=sample_key, y='total', hue=color_by, edgecolor='black', linewidth=0.5, palette=color_palette, dodge=False)
    sns.barplot(data=doublet_counts[doublet_counts['scDblFinder_class'] == 'doublet'], x=sample_key, y='counts', color="#BFBFBF", edgecolor='black', linewidth=0.5)

    plt.ylabel('Nº of cells', size=14, labelpad=10)
    plt.xlabel('')
    plt.xticks(rotation=60, size=12, ha='right')
    plt.legend(frameon = False, loc='best')
    plt.legend("")
    plt.grid(False)

    fig.savefig("QC_doublet_count.png", dpi=300, bbox_inches='tight')
    print("Doublets plot generated!")


def main(adata_dir, sample_key, color_by):

    print('\n>>> Loading adata...', flush=True)
    adata = sc.read_h5ad(adata_dir)
    print('\n>>> adata succesfully loaded', flush=True)

    plot_qc(adata, sample_key=sample_key, color_by= color_by, figsize= None, custom_order = None)
    plot_doublets(adata, sample_key, color_by)

    print("\nQC plots succesfully generetad!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot QC results')
    parser.add_argument('--adata_dir', required=True, help='Path to the adata object from QC module')
    parser.add_argument('--sample_key', required=True, help='Metadata field containing sample names')
    parser.add_argument('--color_by', required=True, help='Metadata field containing the categorical key to color QC plots')
    
    args = parser.parse_args()
    
    main(args.adata_dir,
        args.sample_key,
        args.color_by
    )