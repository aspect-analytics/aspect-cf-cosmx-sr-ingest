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
import dask.array as da
import zarr


output_dir = Path('/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/')
IMAGE_PATH_2D = Path('/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A')
COORD_FP = '/Users/thaotran/GitHub/Project/7_Merck_2024/IBD/CosMx_FF_Morphology2D/IBD7757A/IBD7757A_fov_positions_file.csv.gz'

test_img = TiffFile(glob.glob(str(IMAGE_PATH_2D / "*.TIF"))[0])
ttags = test_img.pages[0].tags
# um/px
image_description = ttags["ImageDescription"].value
image_description_dict = json.loads(image_description)
im_pixel_size_nm = image_description_dict.get("ImPixelSize_nm")
im_pixel_size_nm #120.280945
PIXEL_RES = round(im_pixel_size_nm/1000,3)
# len(image_description_dict.get("Channels")) #5
nb_channels = len(image_description_dict.get("Channels"))
channels_name = [ i['UID'] for i in image_description_dict.get("Channels")]

FOV_SIZE = 4256


### 1. Stich 2D images, try first for DAPI #################################################################
coords = pd.read_csv(COORD_FP, compression='gzip')
coords["x_global_px"] = (coords["X_mm"].values * 10**6 / im_pixel_size_nm).astype(np.uint32)
coords["y_global_px"] = (coords["Y_mm"].values * 10**6 / im_pixel_size_nm).astype(np.int32)

# shift coords
shift_coord_x = max(coords['x_global_px']) - coords['x_global_px']
shift_coord_y = coords['y_global_px']-min(coords['y_global_px'])
coords['shift_coord_x'] = shift_coord_x
coords['shift_coord_y'] = shift_coord_y

max_height = round((max(shift_coord_x)+FOV_SIZE)) 
max_width = round(max(shift_coord_y)+ FOV_SIZE)
# stitched_image = np.zeros((nb_channels,max_width,max_height), dtype=np.uint16) 
stitched_image = np.zeros((max_width,max_height), dtype=np.uint16) 

for i in range(0, len(coords)):
    fov = coords.FOV[i]
    print(fov)
    img_fp = IMAGE_PATH_2D / f"20240327_222234_S2_C902_P99_N99_F{str(fov).zfill(5)}.TIF"
    img = imread (img_fp,key=4)

    x_offset, y_offset = shift_coord_x.iloc[i], shift_coord_y.iloc[i]
    new_x_offset = x_offset+ FOV_SIZE 
    new_y_offset = y_offset + FOV_SIZE 

    stitched_image[y_offset:new_y_offset,
                   x_offset:new_x_offset]=img

# viewer = napari.Viewer()
# viewer.add_image(stitched_image, name="stitched_image", contrast_limits=[0,3000], scale = (1,1))

manager = AspectMicroscopyManager(
    stitched_image,
    "IBD7757A_Cosmx_IF",
    pixel_spacing=(PIXEL_RES,PIXEL_RES),
    reader_func=ArrayLikeReader,
    # reader_kwargs={"channel_names": channels_name}
)
manager.write_scene_zarr(0,
                         output_dir=output_dir,
                         tile_size=2048,
                         compressor="blosc")
