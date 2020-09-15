#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 10:45:24 2020

@author: williamlee
"""

"""
My solutions to the CZI-sponsored "Analysis of single cell RNA-seq data
(Python)" workshop.
"""

##### Tabula Muris #####

import os
path = "/Users/williamlee/scRNA-python-workshop-master"
os.chdir(path)

import pandas as pd 

count_df = pd.read_csv("content/data/brain_counts.csv", index_col=0)
print(count_df.head(2)) # quick check of data
count_df.shape # 3401 cells, 23433 genes

metadata_df = pd.read_csv("content/data/brain_metadata.csv", index_col=0)
print(metadata_df.head(2)) # quick check of metadata
metadata_df.shape # 3401 cells, 5 info. columns

# get # of cells in each subtissue
print(pd.value_counts(metadata_df["subtissue"]))

# summarize each metadata column
meta_cols = metadata_df.columns.values
for col in meta_cols:
    print(pd.value_counts(metadata_df[col]))
    
import scanpy as sc

ann_data = sc.AnnData(X = count_df, obs = metadata_df)
print(ann_data) # quick check of AnnData object

# label spike-ins
spike_ins = {}
num_spike_ins = 0 # instantiate counter
for gene in ann_data.var_names:
    if 'ERCC' in gene:
        spike_ins[gene] = True # spike-in detected
        num_spike_ins += 1
    else:
        spike_ins[gene] = False # spike-in not detected
        
ann_data.var['ERCC'] = pd.Series(spike_ins)
print("number of spike-ins:", num_spike_ins)

# save AnnData object in file (for later)
ann_data.write('content/data/brain_adata.h5ad')

##### Quality Control #####

import matplotlib.pyplot as plt

# read AnnData object from file
ann_data = sc.read('content/data/brain_adata.h5ad')

# get qc metrics object (tuple)
qual_cont = sc.pp.calculate_qc_metrics(ann_data, qc_vars=['ERCC'])

cell_qc_df = qual_cont[0] # store cell qc df
print(cell_qc_df.head(2)) # quick check of cell qc df

gene_qc_df = qual_cont[1] # store gene qc df
print(gene_qc_df.head(2)) # quick check of gene qc df

# plot histogram of # cells vs. total reads (counts)
plt.hist(cell_qc_df['total_counts'], bins=1000)
plt.xlabel('total counts')
plt.ylabel('# cells')
plt.axvline(50000, color='red')
plt.xlim(0,1e6)

# plot histogram of # cells vs. total unique genes
plt.hist(cell_qc_df['n_genes_by_counts'], bins=100)
plt.xlabel('# genes')
plt.ylabel('# cells')
plt.axvline(1000, color='red')

# plot histogram of # cells vs. % ERCC spike-ins (counts)
plt.hist(cell_qc_df['pct_counts_ERCC'], bins=1000)
plt.xlabel('% ERCC counts')
plt.ylabel('# cells')
plt.axvline(10, color='red')

# only keep cells with < 10% ERCC spike-ins
low_ERCC_perc = (cell_qc_df['pct_counts_ERCC']<10)
ann_data = ann_data[low_ERCC_perc]

# filter out cells with < 750 unique genes detected
sc.pp.filter_cells(ann_data, min_genes=750)

# plot histogram of # cells vs. log(# genes)
plt.hist(gene_qc_df['n_cells_by_counts'], bins=1000)
plt.xlabel('# cells (expression > 0)')
plt.ylabel('log(# genes)') 
plt.axvline(2, color='red')
plt.yscale('log')

# plot histogram of total reads (counts) vs. log(# genes)
plt.hist(gene_qc_df['total_counts'], bins=1000)
plt.xlabel('total counts')
plt.ylabel('log(# genes)')
plt.yscale('log') 
plt.axvline(10, color='red')

# filter out genes with < 2 cells expressing or < 10 reads total
sc.pp.filter_genes(ann_data, min_cells = 2)
sc.pp.filter_genes(ann_data, min_counts = 10)

# save quality-controlled AnnData object in file (for later)
ann_data.write('content/data/brain_adata_qc.h5ad')

##### Normalization & PCA #####

ann_data = sc.read('content/data/brain_adata_qc.h5ad')

# perform PCA on quality-controlled AnnData object
sc.pp.pca(ann_data)

# plot PCA colored by mouse ID
sc.pl.pca_overview(ann_data, color='mouse.id')

# cpm normalization (per cell)
ann_data_norm = ann_data.copy() # make copy (for comparison later)
ann_data_norm.raw = ann_data_norm # store copy of raw values
sc.pp.normalize_per_cell(ann_data_norm, counts_per_cell_after=1e6) # normalize

# perform PCA on normalized AnnData object
sc.pp.pca(ann_data_norm)

# plot PCA colored by mouse ID (post-normalization)
sc.pl.pca_overview(ann_data_norm, color='mouse.id')

# cpm normalization (total)
ann_data_norm2 = ann_data.copy()
sc.pp.normalize_total(ann_data_norm2, target_sum=1e6,
                      exclude_highly_expressed=True)
sc.pp.pca(ann_data_norm2) # perform PCA
sc.pl.pca_overview(ann_data_norm2, color='mouse.id') # plot PCA

# remove Rn45s before running PCA 
rm_Rn45s = ann_data_norm.var.index != 'Rn45s'
ann_data_rm_Rn45s = ann_data_norm[:, rm_Rn45s]
sc.pp.pca(ann_data_rm_Rn45s) # perform PCA
sc.pl.pca_overview(ann_data_rm_Rn45s, color='mouse.id') # plot PCA

# centering and scaling gene expression values
sc.pp.log1p(ann_data_norm) # take log(1+x) of values
sc.pp.scale(ann_data_norm) # subtract mean exp. value and divide by SD
sc.pp.pca(ann_data_norm) # perform PCA
sc.pl.pca_overview(ann_data_norm, color='plate.barcode') # plot PCA

# save normalized AnnData object in file (for later)
ann_data_norm.write('content/data/brain_adata_norm.h5ad')

##### Dimensionality reduction #####

ann_data = sc.read('content/data/brain_adata_norm.h5ad')

# run tSNE algorithm on AnnData object
sc.tl.tsne(ann_data, perplexity=30, learning_rate=1000, random_state=0)
sc.pl.tsne(ann_data, color='cell_ontology_class') # plot tSNE

# run UMAP algorithm on AnnData object
sc.pp.neighbors(ann_data) # first compute neighbors graph
sc.tl.umap(ann_data, min_dist=0.5, spread=1.0, random_state=1, n_components=2)
sc.pl.umap(ann_data, color='cell_ontology_class') # plot UMAP

# save normalized AnnData object in file (for later)
ann_data.write('content/data/brain_adata_embed.h5ad')

##### Clustering #####

from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

ann_data = sc.read('content/data/brain_adata_embed.h5ad')

umap_coord = ann_data.obsm['X_umap'] # get UMAP coordinates

# perform k-means clustering using UMAP coordinates
kclust = KMeans(n_clusters=4, random_state=0).fit(umap_coord)

# save k-means labels in obs of AnnData object
ann_data.obs['k-means'] = kclust.labels_
ann_data.obs['k-means'] = ann_data.obs['k-means'].astype(str)

sc.pl.umap(ann_data, color='k-means') # plot k-means

# calculate adjusted Rand index
adj_rand_idx = adjusted_rand_score(labels_true = ann_data.obs['cell_ontology_class'],
                                   labels_pred = ann_data.obs['k-means'])

print("Adjusted Rand index:", round(adj_rand_idx,2)) # print result

# try again with 5 clusters; see if adj. Rand idx. increases
kclust = KMeans(n_clusters=5, random_state=0).fit(umap_coord)
ann_data.obs['k-means'] = kclust.labels_
adj_rand_idx = adjusted_rand_score(labels_true = ann_data.obs['cell_ontology_class'],
                                   labels_pred = ann_data.obs['k-means'])
print("Adjusted Rand index:", round(adj_rand_idx,2))

# utilize louvain algorithm, compare adj. Rand idx. 
sc.tl.louvain(ann_data)
sc.pl.umap(ann_data, color="louvain") # plot louvain
adj_rand_idx = adjusted_rand_score(labels_true = ann_data.obs['cell_ontology_class'],
                                   labels_pred = ann_data.obs['louvain'])
print("Adjusted Rand index:", round(adj_rand_idx,2))

# try again with resolution=0.1; see if adj. Rand idx. increases
sc.tl.louvain(ann_data, resolution=0.1)
adj_rand_idx = adjusted_rand_score(labels_true = ann_data.obs['cell_ontology_class'],
                                   labels_pred = ann_data.obs['louvain'])
print("Adjusted Rand index:", round(adj_rand_idx,2))

# additional exercise
cerebellum = ann_data[ann_data.obs['subtissue'] == 'Cerebellum']
sc.pp.neighbors(cerebellum)
sc.tl.umap(cerebellum)
sc.tl.louvain(cerebellum, resolution=0.1)
sc.pl.umap(cerebellum, color="louvain") # plot louvain
adj_rand_idx = adjusted_rand_score(labels_true = cerebellum.obs['cell_ontology_class'],
                                   labels_pred = cerebellum.obs['louvain'])
print("Adjusted Rand index:", round(adj_rand_idx,2))

# save AnnData object with clusters in file (for later)
ann_data.write('content/data/brain_adata_clusters.h5ad')

##### Differential expression #####

adata = sc.read('content/data/brain_adata_clusters.h5ad')

raw_df = pd.DataFrame(data=adata.raw.X, index=adata.raw.obs_names, columns=adata.raw.var_names)

astr_mark = 'Gja1' # gene of interest here (astrocyte marker)

# here, we are interested in cluster 2
clust2 = raw_df[adata.obs['louvain'] == '2']
rm_clust2 = raw_df[adata.obs['louvain'] != '2']

clust2_markExpr = clust2[astr_mark]
plt.hist(clust2_markExpr.values, bins=100, color='blue', alpha=0.7, label='Cluster 2')
rm_clust2_markExpr = rm_clust2[astr_mark]
plt.hist(rm_clust2_markExpr.values, bins=100, color='red', alpha=0.7, label='Not cluster 2')

plt.ylim(0,100)
plt.xlabel('%s expression'%astr_mark) 
plt.ylabel('# cells')
plt.legend()

from scipy.stats import ttest_ind
ttest_results = ttest_ind(clust2_markExpr, rm_clust2_markExpr,
                          equal_var=False, nan_policy='omit')
print(ttest_results)

# additional exercise
neur = raw_df[adata.obs['cell_ontology_class'] == 'neuron']
endo = raw_df[adata.obs['cell_ontology_class'] == 'endothelial cell']
neur_markExpr = neur['Pecam1']
endo_markExpr = endo['Pecam1']
ttest_results = ttest_ind(neur_markExpr, endo_markExpr,
                          equal_var=False, nan_policy='omit')
print(ttest_results)

# group by cell ontology class
sc.tl.rank_genes_groups(adata, groupby='cell_ontology_class', use_raw=True, 
                        method='t-test_overestim_var', n_genes=10)
sc.pl.rank_genes_groups_tracksplot(adata, groupby='cell_ontology_class')

# group by louvain cluster
sc.tl.rank_genes_groups(adata, groupby='louvain', use_raw=True, 
                        method='t-test_overestim_var', n_genes=10)
sc.pl.rank_genes_groups_tracksplot(adata, groupby='louvain')

# generate heatmap
mark_gene = {
'astrocyte': ['Aldh1l1', 'Slc1a3', 'Aqp4'], 
'oligodendrocyte': ['Mog','Mag'],
'oligodendrocyte precursor cell': ['Pdgfra','Susd5','Cspg4'],
'endothelial cell': ['Pecam1','Cldn5','Slco1c1','Ocln'],
'Bergmann glial cell': ['Gdf10','Vim','Nbl1','A2m'],
'excitatory neuron': ['Slc17a7','Neurod6','Mab21l1'],
'inhibitory neuron': ['Gad1','Reln','Calb1'],
'brain pericyte': ['Des','Mcam','Pdgfrb']
}
sc.pl.matrixplot(adata, mark_gene, groupby='louvain', use_raw=False)
sc.pl.matrixplot(adata, mark_gene, groupby='cell_ontology_class', use_raw=False)