import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import glob
from tifffile import TiffFile
from pathlib import Path
import json
import pandas as pd
import numpy as np
from tifffile import imread,imwrite
from aspect_microscopy_io.manager import AspectMicroscopyManager
from aspect_microscopy_io.reader.array_like_reader import ArrayLikeReader
import napari
import logging
from aspect_merck_toolbox.registration.utils import load_wsireg_omezarr, load_ome_zarr_spacing
import anndata
from anndata import AnnData
import scipy.sparse as sp

output_dir = Path('/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/')
IMAGE_PATH_2D = Path('/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A')
COORD_FP = '/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A/IBD7757A_fov_positions_file.csv.gz'
METADATA_FP = '/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A/IBD7757A_metadata_file.csv.gz'
EXP_MAT_FP = '/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A/IBD7757A_exprMat_file.csv.gz'
FOV_SIZE = 4256
dapi_img_fp = '/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_Cosmx_IF-scene-000.ome.zarr'
PIXEL_RES = load_ome_zarr_spacing(dapi_img_fp)[0]


metadata = pd.read_csv(METADATA_FP, compression='gzip')
exp_mat = pd.read_csv(EXP_MAT_FP, compression='gzip')

if (metadata['cell'] == exp_mat['cell']).all():
    filter_exp_mat = exp_mat.iloc[:, 3:-1]
    filter_exp_mat = filter_exp_mat.drop(filter_exp_mat.filter(regex="^SystemControl").columns, axis=1)
    filter_exp_mat = filter_exp_mat.drop(filter_exp_mat.filter(regex="^Negative").columns, axis=1)
    filter_exp_mat #364492 x 6175 
else:
    print('cell_id not match')


filter_exp_mat.index = (filter_exp_mat.index +1).astype(str)
adata = AnnData(
    X=filter_exp_mat.astype(np.float32),
    obs=metadata,
    var=pd.DataFrame(index=filter_exp_mat.columns),
    # obsm={"spatial": spatial_xy}, ### add this later
    )

adata.X = sp.csr_matrix(adata.X)
adata.write(output_dir/'IBD7757A_Cosmx_raw_anndata.h5ad')

