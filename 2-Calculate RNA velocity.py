import numba #0.52.0
import numpy as np #1.19
import scvelo as scv
import os
import pandas as pd
import sys
from scipy.sparse import csr_matrix
#Input loom file####
loom_data = scv.read("N/velocyto/N.loom")#the file is generated by velocyto (version 0.17.17) with function “run10x”
sample_one = loom_data[np.isin(loom_data.obs.index, sample_obs)]#sample_obs is a list including cell barcodes of cells in specific path
sample_one.obsm['X_umap'] = umap_ordered.values#umap_ordered.values is the coordinate in umap space
adata = sample_one
adata.var_names_make_unique()
#Calculate RNA velocity####
scv.pp.filter_and_normalize(adata, min_shared_counts=10, n_top_genes=3000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=20)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap')
#Dynamical model
scv.pp.filter_and_normalize(adata, min_shared_counts=10, n_top_genes=3000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=20)
scv.tl.recover_dynamics(adata, n_jobs=100)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap')
scv.tl.latent_time(adata)
