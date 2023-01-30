import gseapy as gp
from gseapy import Biomart
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.cluster import hierarchy
from scipy.stats import rankdata
import seaborn as sns
from skimage import filters

#Convenience functions for plotting gene expression matrices and dot plots

def plot_marker_matrix(
    adata, group_key, genes, order=None, gene_order=None,
    cluster_genes=False, cluster_groups=False, otsu_threshold=False,
    annot=False, figsize=None):
  expression_matrix = sc.get.obs_df(adata, [group_key] + list(genes), use_raw=True)
  mean_expression_matrix = expression_matrix.groupby(group_key).mean()

  if otsu_threshold:
    otsu_threshold = mean_expression_matrix.apply(
        filters.threshold_otsu, axis=0)
    plot_matrix = (mean_expression_matrix >= otsu_threshold).astype(int)
    max_expr = expression_matrix.set_index(group_key).max()
    if (max_expr == 0).any():
      print('Dropping genes with 0 expression')
      plot_matrix = plot_matrix.loc[:, max_expr > 0]
    metric = 'hamming'
  else:
    p99 = expression_matrix.set_index(group_key).quantile(0.99)
    sparsity_mask = (p99 == 0)
    p99[sparsity_mask] = expression_matrix.set_index(group_key).max()[sparsity_mask]
    p1 = expression_matrix.set_index(group_key).quantile(0.01)
    plot_matrix = mean_expression_matrix.clip(p1, p99, axis=1)
    if (p99 == 0).any():
      print('Dropping genes with 0 expression')
      plot_matrix = plot_matrix.loc[:, p99 > 0]
    metric = 'correlation'

  if order is not None:
    try:
      plot_matrix = plot_matrix.loc[order]
      mean_expression_matrix = mean_expression_matrix.loc[order]
    except:
      plot_matrix = plot_matrix.iloc[order]
      mean_expression_matrix = mean_expression_matrix.iloc[order]
    cluster_groups = False

  if gene_order is not None:
    plot_matrix = plot_matrix.loc[:, gene_order]
    mean_expression_matrix = mean_expression_matrix.loc[:, gene_order]
    # cluster_genes = False

  if annot:
    annot = mean_expression_matrix
  else:
    annot = None
    
  dendrogram_size = 0.6

  if figsize is None:
    height, width = plot_matrix.shape
    height *= 0.4
    width *= 0.4
    figsize = (width + dendrogram_size, height + dendrogram_size)
  else:
    width, height = figsize

  g = sns.clustermap(plot_matrix, cmap='viridis',
                     cbar_pos=None,
                     dendrogram_ratio=(dendrogram_size / width, dendrogram_size / height),
                     metric=metric,
                     figsize=figsize,
                     yticklabels=True,
                     col_cluster=cluster_genes,
                     row_cluster=cluster_groups,
                     annot=annot,
                     standard_scale=1,
                     vmin=0,
                     fmt='0.1f')
  return g
  

def plot_diffexp_genes(adata, group, num_genes=20):
  diff_exp_df = sc.get.rank_genes_groups_df(adata, group=group)
  overexpressed_genes = diff_exp_df['names'].values[:num_genes]
  underexpressed_genes = diff_exp_df['names'].values[-num_genes:]
  
  plt.figure(figsize=(6, 2))
  sc.pl.rank_genes_groups_violin(
    adata,
    strip=False,
    groups=[group],
    gene_names=overexpressed_genes)

  plt.figure(figsize=(6, 2))
  sc.pl.rank_genes_groups_violin(
    adata,
    strip=False,    
    groups=[group],
    gene_names=underexpressed_genes)

def plot_diffexp_plots(adata, groupby, key, geneset_to_plot, genes_to_plot,
                       plot_extreme_genes=False,
                       num_genes=20, rankgenes_figsize=(6, 4), alpha=0.05):
  plots = {}
  sc.tl.rank_genes_groups(adata, groupby=groupby)
  with plt.rc_context({"figure.figsize": rankgenes_figsize, "axes.grid.axis": 'y'}):
    sc.pl.rank_genes_groups(
        adata, groupby=groupby, sharey=False, ncols=6, n_genes=num_genes,
        fontsize=12)
 
  if plot_extreme_genes:
    plot_diffexp_genes(adata, key)

  diff_exp_df = sc.get.rank_genes_groups_df(adata, group=key)
  rnk = diff_exp_df[['names', 'scores']]

  pre_res = gp.prerank(
      rnk=rnk, gene_sets={geneset_to_plot: genes_to_plot},
      min_size=1, processes=4, permutation_num=10000, seed=0, max_size=10000)

  terms = pre_res.res2d.Term.values
  gsea_plot = gp.plot.gseaplot(rank_metric=pre_res.ranking, term=terms[0],
                               **pre_res.results[terms[0]])
    
  mask = diff_exp_df['names'].isin(genes_to_plot)
  mask &= diff_exp_df['pvals_adj'] < alpha
  return diff_exp_df, gsea_plot

def get_diffexp_matrix(adata, groupby, key, alpha=0.05):
  sc.tl.rank_genes_groups(adata, groupby=groupby)
  diff_exp_df = sc.get.rank_genes_groups_df(adata, group=key)
  mask = diff_exp_df['pvals_adj'] < alpha
  return diff_exp_df[mask]

def plot_ridge_plot(df, group_key='immune_subtype', value_key='log_fc',
                    palette=None, sort=True, xlim=None):
  sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

  # Initialize the FacetGrid object
  if palette is None:
    palette = sns.cubehelix_palette(10, rot=-.25, light=.7)

  if sort:
    group_order = df.groupby(group_key)[value_key].median().sort_values().index
  else:
    group_order = None
  g = sns.FacetGrid(df, row=group_key, hue=group_key, aspect=8,
                    sharey=False, height=.5, palette=palette,
                    row_order=group_order)

  # Draw the densities in a few steps
  g.map(sns.kdeplot, value_key,
        bw_adjust=.5, clip_on=False,
        fill=True, alpha=1, linewidth=1.5)
  g.map(sns.kdeplot, value_key, clip_on=False, color="k", lw=2, bw_adjust=.5)

  # passing color=None to refline() uses the hue mapping
  # g.map(plt.axhline, y=0, linewidth=2, linestyle="-", color=None, clip_on=False)  
  # g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

  if xlim is not None:
    g.set(xlim=xlim)

  # Define and use a simple function to label the plot in axes coordinates
  def label(x, color, label):
      ax = plt.gca()
      txt = ax.text(0, .2, label, fontweight='bold', color=color, fontsize=12,
                    ha="right", va="center", transform=ax.transAxes)
      # txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='#333')])

  g.map(label, group_key)

  # Set the subplots to overlap
  g.fig.subplots_adjust(hspace=-.15)

  # Remove axes details that don't play well with overlap
  g.set_titles("")
  g.set(yticks=[], ylabel="")
  g.despine(bottom=True, left=True)
  return g


def plot_clustered_dot_plot(adata, genes, key, largest_dot=300, figsize=None,
                            categories_order=None, filter_max_expression=True,
                            min_prop_expresed=0.1, keep_markers=None):
  expression_matrix = sc.get.obs_df(adata, list(genes), use_raw=True)  
  if filter_max_expression:
    prop_expressed = expression_matrix.groupby(adata.obs[key].values).apply(lambda x: (x > 0).mean())
    mask = (prop_expressed.max() > min_prop_expresed).values
    if keep_markers is not None:
      mask |= prop_expressed.columns.isin(keep_markers)
    expression_matrix = expression_matrix.loc[:, mask]
    genes = np.array(genes)[mask]
    genes_to_plot = genes[order_genes(expression_matrix)]
  else:
    mask = expression_matrix.max() > 0
    expression_matrix = expression_matrix.loc[:, mask]
    genes_to_plot = np.array(genes)[order_genes(expression_matrix)]
  

  if figsize is None:
    width = len(genes_to_plot)
    height = adata.obs[key].nunique()
    height *= 0.4
    width *= 0.4
    max_dim = np.max([width, height])
    scaling = 1
    if height < 3:
      scaling = 3 / height
    if max_dim > 40:
      scaling = np.max([scaling, 40 / max_dim])
    figsize = (width * scaling, height * scaling)

  g = (sc.pl.DotPlot(adata, genes_to_plot, key, cmap='Reds', figsize=figsize, categories_order=categories_order)
            .style(largest_dot=largest_dot))
  return g


def leiden_to_markers(marker_definitions, thresholded_mean_expression_matrix):
  indicators = []
  n_clusters = thresholded_mean_expression_matrix.shape[0]
  subset_order = []
  for subset, markers in marker_definitions.items():
    subset_indicator = np.ones(n_clusters).astype(bool)
    for gene in markers['+']:
      subset_indicator &= thresholded_mean_expression_matrix[gene].values.astype(bool)
    if '-' in markers:
      for gene in markers['-']:
        subset_indicator &= ~(thresholded_mean_expression_matrix[gene].values.astype(bool))
    indicators.append(subset_indicator)
    subset_order.append(subset)
  indicators = np.stack(indicators, axis=1)  # Leiden clusters by subset
  indicators = pd.DataFrame(
      indicators, index=thresholded_mean_expression_matrix.index,
      columns=subset_order)
  def _decode_rows(row):
    n_matches = np.sum(row > 0)
    if n_matches > 1:
      clusters = ', '.join(row[row > 0].index.values)
      print(f'Leiden cluster {row.name} is assigned to more than 1 cluster ({clusters})')
      raise ValueError()
    if n_matches == 0:
      return None
    return row.index.values[row > 0][0]

  return indicators.apply(_decode_rows, axis=1)


def threshold_markers(mean_expression_matrix, thresholds):
  binarized_markers = {}
  for marker, threshold in thresholds.items():
    binarized_markers[marker] = mean_expression_matrix[marker] > threshold
  return pd.DataFrame(binarized_markers)


def order_genes(gene_expression_matrix, metric='correlation'):
  z = hierarchy.linkage(gene_expression_matrix.T, metric=metric,
                        optimal_ordering=True)
  return hierarchy.leaves_list(z)


def simes_test(df, covariate):
  def _simes(pvals):
    sorted = pvals.sort_values()
    n = len(pvals)
    ranks = np.arange(1, n + 1)
    simes_pval = np.min(sorted * n / ranks)
    return simes_pval    
  return _simes(df[(covariate, "PValue")])    


def fdr(p_vals):
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    return fdr

bm = Biomart(host='uswest.ensembl.org')
h2m = bm.query_simple(dataset='hsapiens_gene_ensembl',
                      attributes=['ensembl_gene_id','external_gene_name',
                                  'mmusculus_homolog_ensembl_gene',
                                  'mmusculus_homolog_associated_gene_name'])
HUMAN_TO_MOUSE_ORTHOLOG_MAP = (h2m.rename(columns={'external_gene_name': 'Gene name', 'mmusculus_homolog_associated_gene_name': 'Mouse gene name'})
                                  .drop_duplicates(subset=['Gene name', 'Mouse gene name']).set_index('Gene name').dropna())

def map_human_to_mouse_orthologs(human_genes):
  missing_human_genes = [g for g in human_genes if g not in HUMAN_TO_MOUSE_ORTHOLOG_MAP.index]
  guesses = [f'{g[0].upper()}{g[1:].lower()}' for g in missing_human_genes]
  subset = [g for g in human_genes if g in HUMAN_TO_MOUSE_ORTHOLOG_MAP.index]
  mouse_genes = HUMAN_TO_MOUSE_ORTHOLOG_MAP.loc[subset]['Mouse gene name']
  return np.hstack([mouse_genes.values, guesses])
