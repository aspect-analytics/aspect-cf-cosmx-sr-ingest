import pandas as pd
import anndata
from anndata import AnnData
import numpy as np
import scanpy as sc
import umap
import h5py
from scipy.spatial.distance import pdist
import numpy as np
from sklearn.neighbors import BallTree
from scipy.sparse import lil_matrix
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix

### Function to calculate boolean adjacency matrix within a threshold
def compute_adjacency_matrix(coordinates, threshold):
    # Build a KDTree for efficient neighborhood queries
    tree = cKDTree(coordinates)
    
    # Query all points within the given threshold
    adj_list = tree.query_pairs(r=threshold)
    
    # Extract row, col indices from the adjacency list
    rows, cols = zip(*adj_list)
    
    # Create the adjacency matrix as a sparse boolean matrix
    adj_matrix = coo_matrix((np.ones(len(rows), dtype=bool), (rows, cols)),
                            shape=(coordinates.shape[0], coordinates.shape[0]))
    
    # Convert the matrix to symmetric by adding its transpose (since it's undirected)
    adj_matrix = adj_matrix + adj_matrix.T
    
    return adj_matrix


### Test with small data
adata_fp = "/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_Cosmx_raw_anndata.h5ad"
raw_adata = anndata.read_h5ad(adata_fp)
coordinates = raw_adata.obsm['spatial'] #364492 × 2
subset_coords = coordinates[:50000]

distances_1 = pdist(subset_coords, metric='euclidean')
boolean_distances_1 = (distances_1 < 30).astype(np.uint8)
adjacency_matrix_1 =  squareform(boolean_distances_1)
adj_matrix_csr1 = csr_matrix(adjacency_matrix_1)

adj_matrix_csr2 = compute_adjacency_matrix(subset_coords, threshold=30)
adjacency_matrix_2 = adj_matrix_csr2.toarray()

# np.array_equal(adjacency_matrix_1,adjacency_matrix_2)
(adj_matrix_csr1 != adj_matrix_csr2).nnz 
scipy.sparse.save_npz(f'/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/Merck_anndata/{sample}/{sample}_{threshold}.npz', adj_matrix_csr2)


### Run for all samples:
import glob
from pathlib import Path
import os
import scipy
adatas_fp = "/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/Merck_anndata"


for fp in glob.glob(adatas_fp + '/*.h5ad'):
    adata = anndata.read_h5ad(fp)
    sample = Path(fp).stem.split('.')[0]
    print(sample)
    coordinates = adata.obsm['spatial']*0.12 # if coordinate is pixel, convert to micron
    for threshold in [15, 20, 25, 30, 35, 40,50]:
        print(threshold)
        adj_matrix_csr = compute_adjacency_matrix(coordinates, threshold=threshold)
        file_path = f'/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/Merck_anndata/{sample}/{sample}_{threshold}.npz'
        folder = os.path.dirname(file_path)
        os.makedirs(folder, exist_ok=True) 
        scipy.sparse.save_npz(file_path, adj_matrix_csr)









