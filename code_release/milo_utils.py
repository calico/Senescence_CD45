# Convenience functions for Milo differential analysis based plots
import scanpy as sc
import matplotlib.pyplot as plt
import milopy.core as milo
import milopy.plot as milopl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
import re
import rpy2
from rpy2.robjects.packages import PackageNotInstalledError,importr
import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import STAP
import seaborn as sns

colors = sns.color_palette('colorblind', as_cmap=True) + sns.color_palette('dark', as_cmap=True) + sns.color_palette('bright', as_cmap=True)

# MOSTLY COPIED FROM MILOPY WITH MODIFICATIONS

def _try_import_bioc_library(name):
    try:
        _r_lib = importr(name)
        return _r_lib
    except PackageNotInstalledError:
        raise RuntimeError(
            f"Install Bioconductor library `{name!r}` first as `BiocManager::install({name!r}).`"
        )

def set_reference_level(df, col, reference_level):
    '''R assumes the first categorical level is the reference, this function
       re-orders a categorical variable so that the desired reference is the
       first categorical level.
    '''
    new_categories = [reference_level] + [x for x in df[col].cat.categories if x != reference_level]
    df.loc[:, col] = df[col].cat.reorder_categories(new_categories)
    return df


def _graph_spatialFDR(res, kth_distance):
    '''
    FDR correction weighted on inverse of connectivity of neighbourhoods.
    The distance to the k-th nearest neighbor is used as a measure of connectivity.
    '''  
    # use 1/connectivity as the weighting for the weighted BH adjustment from Cydar
    w = 1 / kth_distance
    w = w.loc[res.index]
    w[np.isinf(w)] = 0

    ## Computing a density-weighted q-value.
    pvalues = res["PValue"]
    o = pvalues.argsort()
    pvalues = pvalues[o]
    w = w[o]

    adjp = np.zeros(shape=len(pvalues))
    adjp[o] = (sum(w) * pvalues / np.cumsum(w))[::-1].cummin()[::-1]
    adjp = np.clip(np.array(adjp), 0, 1)

    res.loc[:, 'SpatialFDR'] = adjp
    return res

## Set up rpy2 to run edgeR

def DA_nhoods(adata, design, reference_levels=None, subset_samples=None,
              annotation_cols=None, embedding_key='X_umap'):
  rpy2.robjects.numpy2ri.activate()
  rpy2.robjects.pandas2ri.activate()
  edgeR = _try_import_bioc_library("edgeR")
  limma = _try_import_bioc_library("limma")
  stats = importr("stats")
  base = importr("base")

  nhood_adata = adata.uns["nhood_adata"]
  covariates = [x.strip(" ") for x in set(re.split('\\+|\\*', design.lstrip("~ ")))]

  ## Add covariates used for testing to nhood_adata.var
  sample_col = nhood_adata.uns["sample_col"]
  try:
      nhoods_var = adata.obs[covariates + [sample_col]].drop_duplicates()
  except KeyError:
      missing_cov = [x for x in covariates if x not in nhood_adata.var.columns]
      raise KeyError(
          'Covariates {c} are not columns in adata.obs'.format(c=" ".join(missing_cov))
      )
  ## N.B. This might need some type adjustment!!
  nhoods_var = nhoods_var[covariates + [sample_col]]
  nhoods_var.index = nhoods_var[sample_col].astype("str")

  try:
      assert nhoods_var.loc[nhood_adata.var_names].shape[0] == len(nhood_adata.var_names)
  except:
      raise ValueError("Covariates cannot be unambiguously assigned to each sample -- each sample value should match a single covariate value")
  nhood_adata.var = nhoods_var.loc[nhood_adata.var_names]

  ## Get design dataframe
  try:
      design_df = nhood_adata.var[covariates].copy()
  except KeyError:
      missing_cov = [x for x in covariates if x not in nhood_adata.var.columns]
      raise KeyError(
          'Covariates {c} are not columns in adata.uns["nhood_adata"].var'.format(c=" ".join(missing_cov))
      )

  # Set reference levels for categorical covariates (aka factors)
  if reference_levels is not None:
    for col, reference_level in reference_levels.items():
      design_df = set_reference_level(design_df, col, reference_level)

  ## Get count matrix
  count_mat = nhood_adata.X.toarray()
  lib_size = count_mat.sum(0)

  ## Filter out samples with zero counts
  keep_smp = lib_size > 0

  ## Subset samples
  if subset_samples is not None:
      keep_smp = keep_smp & nhood_adata.var_names.isin(subset_samples)

  ## Filter out nhoods with zero counts
  ## (they can appear after sample filtering)
  keep_nhoods = count_mat[:,keep_smp].sum(1) > 0

  model = stats.model_matrix(object=stats.formula(design), data=design_df[keep_smp])

  dge = edgeR.DGEList(counts=count_mat[keep_nhoods,:][:,keep_smp])
  dge = edgeR.calcNormFactors(dge, method="TMM")
  dge = edgeR.estimateDisp(dge, model)
  fit = edgeR.glmQLFit(dge, model, robust=True)

  # Get model col names
  r_str = '''
  get_model_cols <- function(design_df, design){
      m = model.matrix(object=formula(design), data=design_df)
      return(colnames(m))
  }
  '''
  get_model_cols = STAP(r_str, "get_model_cols")
  model_mat_cols = get_model_cols.get_model_cols(design_df, design)
  n_coef = len(model_mat_cols)
  assert n_coef == model.shape[1]

  results = []

  # MODIFIED: Loop through to get results for all coefs instead of just
  # taking the last coef
  for i, col_name in enumerate(model_mat_cols):
    if col_name == '(Intercept)':
      continue
    res = base.as_data_frame(edgeR.topTags(edgeR.glmQLFTest(fit, coef=i+1), sort_by='none', n=np.inf))
    # Different versions of edgeR / rpy2 seem to return res in different formats.
    res = pd.DataFrame(res)
    if 'PValue' not in res.columns:
        res = res.T
        res.columns = ['logFC', 'logCPM', 'F', 'PValue', 'FDR']
    res.index = nhood_adata.obs_names[keep_nhoods]
    res = _graph_spatialFDR(res, nhood_adata.obs.loc[keep_nhoods, 'kth_distance'])
    res.columns = [(col_name, feature) for feature in res.columns]
    results.append(res)

  results_df = pd.concat(results, axis=1)
  results_df.columns = pd.MultiIndex.from_tuples(
      results_df.columns, names=['coef', ''])
  
  # MODIFIED: Add some basic metadata for each nhood for convenient plotting downstream
  results_df.loc[:, 'index_cell'] = nhood_adata.obs['index_cell'].values
  results_df.loc[:, 'kth_distance'] = nhood_adata.obs['kth_distance'].values
  embedding_coords = pd.DataFrame(adata.obsm[embedding_key],
                                  index=adata.obs.index,
                                  columns=['x', 'y'])
  embedding_coords = embedding_coords.loc[results_df['index_cell'].values]
  results_df.loc[:, 'x'] = embedding_coords['x'].values
  results_df.loc[:, 'y'] = embedding_coords['y'].values

  # Assign nhoods to annotations
  if annotation_cols is not None:
    for annotation_col in annotation_cols:
      annotations = adata.obs[annotation_col]
      assignments = []
      for i, nhood_indicator in enumerate(adata.obsm['nhoods'].T):
        indicators = np.asarray(nhood_indicator.todense()).flatten().astype(bool)
        assignments.append(annotations[indicators].value_counts().idxmax())
      assignments = np.array(assignments)
      results_df.loc[:, annotation_col] = assignments
  return results_df


def plot_nhood_da_results(adata, results_df, covariate, annotation_col,
                          significance=0.1, figsize=(12, 12), palette=None):
  plt.figure(figsize=figsize)  
  logFC = results_df[covariate]['logFC']
  spatialFDR = results_df[covariate]['SpatialFDR']
  mask = spatialFDR < significance

  # Plot t-SNE
  ax = plt.subplot(1, 3, 1)
  ax = sc.pl.umap(
      adata,
      color='p53_condition',
      frameon=False,
      palette=colors,
      wspace=0.1,
      legend_fontsize='small',
      ax=ax,
      show=False)
    
  scp = plt.scatter(results_df[mask]['x'],
                   results_df[mask]['y'],
                   c=logFC[mask],
                   s=150*np.abs(logFC[mask]),
                   lw=1,
                   edgecolors='k',
                   alpha=0.75,
                   cmap='RdBu_r')

  axins = inset_axes(ax,
                     width="100%",  
                     height="5%",
                     loc='lower center',
                     borderpad=-5)
  plt.colorbar(scp, cax=axins, orientation="horizontal")

  ax = plt.subplot(1, 3, 2)
  ax = sc.pl.umap(
      adata,
      color=[annotation_col],
      frameon=False,
      palette=palette,
      wspace=0.1,
      legend_fontsize='small',
      ax=ax,
      show=False)
  plt.axis('off')
  plt.grid(b=None)
  sns.despine(left=True, bottom=True)

  ax = plt.subplot(1, 3, 3)
  sorted_annos = results_df.groupby(annotation_col).\
      median().\
      sort_values((covariate, "logFC")).index

  anno_df = pd.concat([results_df[annotation_col], results_df[covariate]], axis=1)
  anno_df.loc[:, 'is_signif'] = mask
  sns.violinplot(data=anno_df, y=annotation_col, x="logFC", order=sorted_annos, 
                 size=190, 
                 inner=None,
                 orient="h",
                 color='#aaa',
                 linewidth=0,
                 ax=ax,
                 scale="area")

  sns.stripplot(data=anno_df, y=annotation_col, x="logFC", order=sorted_annos,
                size=5,
                hue='is_signif',
                palette={True: 'r', False: '#555'},
                ax=ax,
                orient="h",
                alpha=1.0);
  plt.legend(['not significant', 'significant'], facecolor="white")
  leg = plt.gca().get_legend()
  leg.legendHandles[0].set_color('#555')
  leg.legendHandles[1].set_color('r')
  plt.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--")
  yticklabels = plt.gca().get_yticklabels()
  for label in yticklabels:
    label.set_color(palette[label.get_text()])
  plt.tight_layout()


def add_enrichment_annotation(adata, results_df, covariate, significance=0.2):
  covariate_results = results_df[covariate]
  neighborhood_is_sig = covariate_results['SpatialFDR'] < significance
  neighborhood_is_sig_enriched = neighborhood_is_sig & (covariate_results['logFC'] > 0)
  neighborhood_is_sig_depleted = neighborhood_is_sig & (covariate_results['logFC'] < 0)

  membership_matrix = adata.obsm['nhoods'].toarray()
  cell_in_sig_enriched_neighborhood = membership_matrix[:, neighborhood_is_sig_enriched].sum(axis=1) > 0
  cell_in_sig_depleted_neighborhood = membership_matrix[:, neighborhood_is_sig_depleted].sum(axis=1) > 0
  
  adata.obs['in_sig_enriched_milo_neighborhood'] = pd.Series(
      cell_in_sig_enriched_neighborhood, index=adata.obs.index).astype('category')
  adata.obs['in_sig_depleted_milo_neighborhood'] = pd.Series(
      cell_in_sig_depleted_neighborhood, index=adata.obs.index).astype('category')
  return adata


def run_milo_analysis(adata, sample_col, design,
                      reference_levels=None,
                      n_neighbors=40,
                      n_pcs=10,
                      annotation_cols=None):
  ## Run PCA  
  sc.tl.pca(adata, svd_solver='arpack')

  ## Build KNN graph
  sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

  ## Assign cells to neighbourhoods
  milo.make_nhoods(adata)

  ## Count cells from each sample in each nhood
  milo.count_nhoods(adata, sample_col=sample_col)

  ## Test for differential abundance between conditions
  results_df = DA_nhoods(adata,
                         design=design,
                         reference_levels=reference_levels,
                         annotation_cols=annotation_cols)
  return results_df