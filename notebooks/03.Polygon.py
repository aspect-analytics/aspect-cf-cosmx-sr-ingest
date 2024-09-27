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
from aspect_merck_toolbox.visualization.utils import load_omezarr_pyr
import shapely
from shapely.geometry import Polygon
from shapely.affinity import scale
import cv2
import dask.array as da
import anndata
import zarr


output_dir = Path('/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/')
IMAGE_PATH_2D = Path('/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A')
COORD_FP = '/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A/IBD7757A_fov_positions_file.csv.gz'
METADATA_FP = '/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A/IBD7757A_metadata_file.csv.gz'
EXP_MAT_FP = '/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A/IBD7757A_exprMat_file.csv.gz'
CELL_POLYGON_FP = '/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A/IBD7757A-polygons.csv.gz'
TX_FILE = '/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A/IBD7757A_tx_file.csv.gz'
FOV_SIZE = 4256
dapi_img_fp = '/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_Cosmx_IF-scene-000.ome.zarr'
PIXEL_RES = load_ome_zarr_spacing(dapi_img_fp)[0]
dapi_img = load_omezarr_pyr(dapi_img_fp, move_axis=False)[0]

def get_im_pixel_size_nm (img_path):
    img = TiffFile(glob.glob(str(img_path / "*.TIF"))[0])
    ttags = img.pages[0].tags
    image_description = ttags["ImageDescription"].value
    image_description_dict = json.loads(image_description)
    im_pixel_size_nm = image_description_dict.get("ImPixelSize_nm")
    return im_pixel_size_nm


im_pixel_size_nm = get_im_pixel_size_nm(IMAGE_PATH_2D)
metadata = pd.read_csv(METADATA_FP, compression='gzip')
cell_polygon = pd.read_csv(CELL_POLYGON_FP, compression='gzip')
coords = pd.read_csv(COORD_FP, compression='gzip')


# Shift FOV
coords["x_global_px"] = (coords["X_mm"].values * 10**6 / im_pixel_size_nm).astype(np.uint32)
coords["y_global_px"] = (coords["Y_mm"].values * 10**6 / im_pixel_size_nm).astype(np.int32)
shift_coord_x = max(coords['x_global_px']) - coords['x_global_px']
shift_coord_y = coords['y_global_px']-min(coords['y_global_px'])
coords['shift_coord_x'] = shift_coord_x
coords['shift_coord_y'] = shift_coord_y

# Shift polygon
shift_coord_x = max(cell_polygon['x_global_px']) - cell_polygon['x_global_px']
shift_coord_y = cell_polygon['y_global_px'] - min(cell_polygon['y_global_px'])
cell_polygon['shift_coord_x'] = shift_coord_x.astype(np.int32)
cell_polygon['shift_coord_y'] = shift_coord_y.astype(np.int32)

# Function to flip polygons 
def flip_polygons_shapely(fov_polygons, center, flip_x=True, flip_y=True):
    # flipped_polys = []
    if flip_x:
        xfact = -1.0
    if flip_y:
        yfact = -1.0
    for idx, cell_data in fov_polygons.groupby("cell"):
        # Create a shapely polygon object from the polygon coordinates
        polygon = Polygon(zip(cell_data["shift_coord_y"], cell_data["shift_coord_x"]))
        # Scale/flip the polygon relative to its centroid
        flipped_polygon = scale(polygon, xfact=xfact, yfact=yfact, origin=center)
        # Append the flipped polygon coordinates as a numpy array
        cell_bounds[idx] = np.array(flipped_polygon.exterior.coords)
        # flipped_polys.append(np.array(flipped_polygon.exterior.coords))
    # return flipped_polys


# Flip polygons in all FOV:
len(coords) #238
cell_bounds = {}
for fov in range(1,len(coords)+1):
    print(fov)
    fov_poly_i = cell_polygon[cell_polygon["fov"] == fov]
    fov_coords_i = coords[coords["FOV"] == fov]
    center = (fov_coords_i["shift_coord_y"].values[0] + 1/2*FOV_SIZE,fov_coords_i["shift_coord_x"].values[0]+ 1/2*FOV_SIZE)
    flip_polygons_shapely(fov_poly_i, center, flip_x=True, flip_y=True)

# Filter cells not in anndata
unique_cell_id = metadata['cell'].tolist()
if unique_cell_id != [*cell_bounds]:
    print("cell id not match")
     # Create a new dictionary with keys arranged in the metadata order
    arranged_cell_bounds = {key: cell_bounds[key] for key in unique_cell_id if key in cell_bounds}

# Bitmask
bitmask_shape = np.array([dapi_img.shape[0],
                          dapi_img.shape[1]])

cell_bitmask = np.zeros(bitmask_shape, dtype=np.int32)
idx = 1
for cell, poly in arranged_cell_bounds.items():
    if np.isnan(poly).any():
        print(f"Skipping cell {cell} due to NaN values")
        continue
    poly = np.asarray([poly], dtype=np.int32)
    poly_swapped = poly.copy()
    poly_swapped[:, :, [0, 1]] = poly_swapped[:, :, [1, 0]]
    cell_bitmask = cv2.fillPoly(cell_bitmask, [poly_swapped], idx)
    idx+=1
cell_bitmask = cell_bitmask.astype(np.uint32)
# imwrite('/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_cell_bitmask.tif', cell_bitmask, compression='deflate')

cell_contour_mask = np.zeros(bitmask_shape, dtype=np.int32)
idx = 1
for cell, poly in arranged_cell_bounds.items():
    if np.isnan(poly).any():
        print(f"Skipping cell {cell} due to NaN values")
        continue
    poly = np.asarray([poly], dtype=np.int32)
    poly_swapped = poly.copy()
    poly_swapped[:, :, [0, 1]] = poly_swapped[:, :, [1, 0]]
    cell_contour_mask = cv2.polylines(
        cell_contour_mask,
        [np.squeeze(poly_swapped).astype(np.int32)],
        True,
        idx,
    )
    idx+=1
cell_contour_mask = cell_contour_mask.astype(np.uint32)
# imwrite('/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_cell_contour_mask.tif', cell_contour_mask, compression='deflate')
# combine = np.stack([cell_bitmask,cell_contour_mask], axis=0).astype(np.uint32)
bm_fp = ['/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_cell_bitmask.tif',
            '/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_cell_contour_mask.tif']
combine = da.stack([da.from_zarr(zarr.open(imread(fp,aszarr=True, series=0))) for fp in bm_fp], axis=0)

manager = AspectMicroscopyManager(
    combine,
    "IBD7757A_Cosmx_Bitmask",
    reader_func=ArrayLikeReader,
    pixel_spacing=(PIXEL_RES,PIXEL_RES),
    reader_kwargs={"channel_names": ["Filled Cells", "Cell Contours"]},
)

manager.write_scene_zarr(
    0,
    output_dir=output_dir,
    tile_size=2048,
    compressor="blosc",
    subsample=True
)


### Put polygon coordinate into anndata #################################################################
cell_centroid = {}
for cell, poly in arranged_cell_bounds.items():
    polygon = Polygon(poly)
    poly_centroid = polygon.centroid
    cell_centroid[cell] = (poly_centroid.y, poly_centroid.x)

spatial_coord = np.array([list(cell_centroid[cell]) for cell in cell_centroid.keys()]).astype(np.float32) * PIXEL_RES
spatial_coord

RAW_ADATA_FP = '/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_Cosmx_raw_anndata.h5ad'
adata = anndata.read_h5ad(RAW_ADATA_FP)
if adata.obs['cell'].tolist() == [*cell_centroid]:
    adata.obsm['spatial'] = spatial_coord
else:
    print("cell id bitmask not match anndata")

adata.X = adata.X.astype(np.float32)
adata.X = adata.X.toarray()
adata_zarr_fp = output_dir/'IBD7757A_Cosmx_raw_anndata.zarr'
adata.write_zarr(adata_zarr_fp, chunks=(adata.shape[0], 10))
