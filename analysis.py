#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :analysis.py
@说明        : just do it!
@时间        :2022/03/19 00:26:31
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com
'''


# 注意：该脚本基本由scanpy和scvelo这两个包的官方教程组成，其中参数复杂多样，故不一一规定，按照个人分析需求，尝试不同参数即可

from mpl_toolkits.axes_grid1 import make_axes_locatable
import scvelo as scv
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

in_adata = 'wwjtest_bin_50.h5ad'
adata = sc.read(in_adata)
adata.X = adata.X.astype('float64')

sc.settings.verbosity = 3
scv.settings.verbosity = 3
figdir = 'wwjtest10'
sc.settings.figdir = figdir
scv.settings.figdir = figdir
sc.settings.set_figure_params(
    dpi_save=300, figsize=[12, 12], facecolor='#dcdfe6', format='png')
scv.set_figure_params(dpi_save=300, figsize=[
                      12, 12], facecolor='#dcdfe6', format='png')

x_cor = adata.obsm['X_spatial'][:, 0]
y_cor = adata.obsm['X_spatial'][:, 1]
x_len = x_cor.max()-x_cor.min()
y_len = y_cor.max()-y_cor.min()


# before filter
# draw info
sc.pl.highest_expr_genes(adata, n_top=30, save='_raw')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=True, inplace=True)
sc.pl.violin(adata, ['n_counts', 'n_genes_by_counts', 'total_counts'],
             jitter=0.4, multi_panel=True, save='_cell_raw')
sc.pl.violin(adata.T, ['n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts',
             'total_counts'], jitter=0.4, multi_panel=True, save='_gene_raw')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_raw')

# filter the data
# filtered out cells that have less than 20 genes expresse
sc.pp.filter_cells(adata, min_genes=50)
sc.pp.filter_cells(adata, max_genes=1000)  # !4000
# filtered out genes that are detected in less than 5 cells
sc.pp.filter_genes(adata, min_cells=3)
#sc.pp.filter_genes(adata, max_cells=1500)
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=True, inplace=True)
# draw info
sc.pl.highest_expr_genes(adata, n_top=30, save='_filter')
sc.pl.violin(adata, ['n_counts', 'n_genes_by_counts', 'total_counts'],
             jitter=0.4, multi_panel=True, save='_cell_filter')
sc.pl.violin(adata.T, ['n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts',
             'total_counts'], jitter=0.4, multi_panel=True, save='_gene_filter')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_filter')
# -----> it should be better!!
adata = adata[adata.obs.total_counts < 60000, :]

# normalize the data
# sc.pp.normalize_total(adata)
sc.pp.normalize_total(adata, target_sum=1e6)  # 1000000reads per cell??  CPM?
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
#sc.pp.highly_variable_genes(adata, max_mean=10)
sc.pp.highly_variable_genes(adata, max_mean=10)  # ok
sc.pl.highly_variable_genes(adata, save='')

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
# sc.pp.scale(adata)
sc.pp.scale(adata, max_value=10)

# pca
sc.tl.pca(adata, svd_solver='arpack', n_comps=30, use_highly_variable=True)
sc.pl.pca(adata, color='sample', save='_sample')
sc.pl.pca_variance_ratio(adata, log=True, save='')
pca_df = pd.DataFrame(adata.obsm['X_pca'])[[0, 1, 2]]
pca_df.to_csv('3pc.csv', index=False, header=None)

# Computing the neighborhood graph use pcs
sc.pp.neighbors(adata, n_neighbors=8, n_pcs=10, use_rep='X_pca')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=15, metric='l1')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10, metric='euclidean')  # $$$$
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=10, metric='manhattan')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=10, metric='correlation')
sc.pp.neighbors(adata, n_neighbors=30, metric='euclidean')  # $$$$ook
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=15, metric='euclidean')  # 10

sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10, use_rep='X_diffmap')

sc.tl.leiden(adata, resolution=1, use_weights=True)
sc.tl.louvain(adata, resolution=1, use_weights=True)
sc.tl.umap(adata, min_dist=0.001)  # ok
sc.tl.umap(adata, min_dist=0.0001, spread=10, maxiter=1000)

sc.pl.umap(adata, color=['sample', 'louvain', 'leiden'],
           legend_loc='on data', save='')

# metric : Union[Literal[‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’, ‘manhattan’],
#                Literal[‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘correlation’, ‘dice’,
#                        ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’,
#                        ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’,
#                        ‘sokalsneath’, ‘sqeuclidean’, ‘yule’],

# sc.tl.diffmap(adata)
#sc.pp.neighbors(adata, n_neighbors=15, n_pcs=5,use_rep='X_diffmap')

sc.tl.draw_graph(adata, n_jobs=50)
sc.pl.draw_graph(adata, color=['sample', 'louvain',
                 'leiden'], legend_loc='right margin', save='_raw')


# spacital
scv.pl.scatter(adata, color='leiden', basis='spatial', alpha=1, legend_loc='right margin',
               size=45, figsize=[12, 12*y_len/x_len], save='spatial_leiden')
scv.pl.scatter(adata, color='louvain', basis='spatial', alpha=1, legend_loc='right margin',
               size=45, figsize=[12, 12*y_len/x_len], save='spatial_louvain')

# find marker genes
# what means wilcoxon
sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes')
sc.pl.rank_genes_groups_heatmap(
    adata, n_genes=10, groupby="louvain", save='_marker_genes_heatmap', show_gene_labels=True)
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
     for group in groups for key in ['names', 'pvals']}) \
    .to_csv('marker.genes.tsv', sep='\t', index=False)
sc.pl.rank_genes_groups_violin(adata, n_genes=25, save='_wilcoxon')

# paga
sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, threshold=0.03, save='')
sc.tl.draw_graph(adata, init_pos='paga', n_jobs=50)
sc.pl.draw_graph(adata, color=['sample', 'louvain',
                 'leiden'], legend_loc='on data', save='_paga')
sc.pl.paga(adata, color=['sample', 'louvain',
           'leiden'], save='_louvain_leiden')
sc.pl.paga_compare(adata, threshold=0.03, title='', right_margin=0.2,
                   size=10, edge_width_scale=0.5, legend_fontsize=12, fontsize=12,
                   frameon=False, edges=True, save=True)

# adata.write_h5ad('wwj.scanpy_ok.h5ad')
adata = sc.read('wwj.scanpy_ok.h5ad')
#adata.X = adata.X.astype('float64')

scv.pl.proportions(adata, dpi=600, groupby='louvain',
                   figsize=[12, 12], save='')
scv.pp.filter_and_normalize(adata, n_top_genes=3000)  # n_top_genes ??
scv.pp.filter_and_normalize(adata, log=False)  # log for two layers? ??  ##ok
scv.pp.filter_and_normalize(adata, max_mean=10, layers_normalize=[
                            'spliced', 'unspliced'])

scv.pp.filter_genes(adata, min_shared_counts=3)  # 10
scv.pp.normalize_per_cell(adata, layers=['spliced', 'unspliced'])  # 10
# if n_top_genes is not None:
#     scv.pp.filter_genes_dispersion(adata,n_top_genes=3000)  #10
# if log:
#     scv.pp.log1p(adata)


scv.pp.moments(adata, use_highly_variable=True)

# scv.tl.recover_dynamics(adata,n_jobs=8,fit_connected_states=True)
scv.tl.recover_dynamics(adata, n_top_genes=3000, n_jobs=8)
scv.tl.recover_dynamics(adata, n_jobs=8)  # 10

#scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity(adata, groupby='sample', mode='dynamical')
scv.tl.velocity_graph(adata, n_jobs=8)

scv.pl.velocity_embedding_stream(adata, basis='umap', save='dyn_stream_umap')
scv.pl.velocity_embedding_stream(adata, basis='pca', save='dyn_stream_pca')
scv.pl.velocity_embedding_stream(
    adata, basis='draw_graph_fa', save='dyn_stream_draw_graph_fa')
scv.pl.velocity_embedding_stream(adata, basis='spatial', alpha=1,
                                 size=45, arrow_size=1, min_mass=2, figsize=[12, 12*y_len/x_len],
                                 color='louvain', legend_loc='right margin', save='dyn_stream_spatial_louvain')
scv.pl.velocity_embedding_stream(adata, basis='spatial', alpha=1,
                                 size=45, arrow_size=1, min_mass=2, figsize=[12, 12*y_len/x_len],
                                 color='leiden', legend_loc='right margin', save='dyn_stream_spatial_leiden')
scv.pl.velocity_embedding(adata, basis='umap', arrow_length=3,
                          arrow_size=2, save='dyn_arrow_umap')
scv.pl.velocity_embedding(adata, basis='pca', arrow_length=3,
                          arrow_size=2, save='dyn_arrow_pca')
scv.pl.velocity_embedding(adata, basis='draw_graph_fa', arrow_length=3,
                          arrow_size=2, save='dyn_arrow_draw_graph_fa')

scv.pl.velocity_embedding(adata, basis='spatial', alpha=0.7,
                          arrow_length=3, arrow_size=2, figsize=[12, 12*y_len/x_len],
                          color='louvain', legend_loc='right margin', save='dyn_arrow_spatial_louvain')
scv.pl.velocity_graph(adata, threshold=.1, save='draw_graph')
scv.tl.paga(adata, groups='louvain')
scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, save='paga')
# calculate and plot the latent time
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', basis='spatial',
               size=45, figsize=[12, 12*y_len/x_len], color_map='gnuplot', save='latent_time')

top_genes = adata.var['fit_likelihood'].sort_values(
    ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time',
               col_color='louvain', n_convolve=100, save='gene_300')
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5,
               figsize=[4, 4], save='top_likehood_genes_all')
marker_genes = ['Zm00001eb072880', 'Zm00001eb335880',
                'Zm00001eb320810', 'Zm00001eb046400']
scv.pl.scatter(adata, marker_genes, save='marker_gene_likehood')
scv.pl.scatter(adata, x='latent_time', y=marker_genes,
               save='marker_gene_likehood_latent_time')

scv.pl.velocity(adata, marker_genes, ncols=2, save='_wwjtest', basis='spatial')
scv.pl.velocity(adata, top_genes[:6], ncols=2,
                save='_wwjtest', basis='spatial')

scv.pl.scatter(adata, color='velocity_pseudotime', save='_pseudotime_scatter')
scv.pl.scatter(adata, color='latent_time', save='_latent_time_scatter')

# discover dynamical genes
scv.tl.rank_dynamical_genes(adata, groupby='louvain')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.to_csv('Cluster-specific_top-likelihood_genes.tsv', sep='\t', index=False)
for clu in df.columns:
    scv.pl.scatter(adata, basis=df[clu][:15], ncols=5, figsize=[4, 4],
                   save='top_likehood_genes_cluster{}'.format(clu))

# save the result in h5ad format
adata.write_h5ad('scvelo_wwjok.h5ad')
marker_genes = adata.var['fit_likelihood'].sort_values(
    ascending=False).index[:20]

sc.pl.dotplot(adata, marker_genes, groupby='louvain',
              save='marker_gene_expre_louvain')
