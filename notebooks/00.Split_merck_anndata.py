import pandas as pd
import anndata
from anndata import AnnData
import numpy as np
import scanpy as sc
import h5py
from scipy.sparse import csr_matrix


adata_fp = '/home/ubuntu/IBD/CosMx_anndata/0903_Dipart_subtyped_nn.h5ad'
celltype_fp = '/home/ubuntu/IBD/CosMx_anndata/0903_Dipart_subtyped_nn.csv'

celltype_df = pd.read_csv(celltype_fp)
celltype_df['index'] # {fov}_{cell_ID}+{Sample_ID} (Sample_ID is in IBD_Commercial_Pilot_CosMx_Metadata(in).csv)
sample_id = [i.split('-')[1] for i in celltype_df['index'].astype(str)]
celltype_df['sample_id'] = sample_id

sample_list = celltype_df['sample_id'].unique().tolist()

for sample_id in sample_list[1:]:
    print(sample_id)
    # 0. Obs + Celltyping
    subset_obs = celltype_df[celltype_df['sample_id'] == sample_id]
    # 1. X_sparse
    subset_indices = celltype_df[celltype_df['sample_id'] == sample_id].index.values  # Get the row indices for this sample
    with h5py.File(adata_fp, 'r') as f:
        X = f['X']  
        X_subset = X[subset_indices, :]  
    X_subset = X_subset.astype('float32')
    X_sparse = csr_matrix(X_subset)
    # 2. Var
    with h5py.File(adata_fp, 'r') as f:
        var = f['var']
        var_df = pd.DataFrame({key: var[key][()] for key in var.keys()})
    var_df
    # 3. Raw
    # with h5py.File(adata_fp, 'r') as f:
    #     raw_subset = f['layers']['raw'][subset_indices, :]  
    # raw_subset= raw_subset.astype('float32')
    # 4.Obsm
    subset_obsm = {}
    obsm_list = ['X_cellcharter', 'X_pca', 'X_pca_original', 'X_scVI', 'X_umap', 'spatial']
    with h5py.File(adata_fp, 'r') as f:
        obsm = f['obsm']
        for key in obsm_list:
            print(key)
            data = pd.DataFrame(obsm.get(key))
            subset_data = data.loc[subset_indices, :]
            subset_obsm[key] = subset_data
    ### Create new anndata object
    adata = AnnData(X=X_sparse, var=var_df, obs=subset_obs)
    # adata.layers['raw'] = raw_subset
    for key in subset_obsm.keys():
        adata.obsm[key] = subset_obsm[key].values.astype('float32')

    adata.obs['predicted_doublet_no_sim'] = adata.obs['predicted_doublet_no_sim'].astype('boolean')
    adata.var.index = adata.var['_index'].astype(str).values
    adata.var = adata.var.drop(columns=['_index'])

    adata.write_h5ad(f'/home/ubuntu/IBD/CosMx_anndata/{sample_id}.h5ad')




