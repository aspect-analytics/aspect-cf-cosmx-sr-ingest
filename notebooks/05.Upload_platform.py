from pathlib import Path
import pandas as pd
from aspect_sdk.platform_sdk.client import SDKClient

AA_DEV = "https://aa-dev.platform.aspect-analytics.com/"
sdk_client = SDKClient.create(url=AA_DEV)

# Bitmask
bitmas_up = sdk_client.asset_storage_client.upload_new_asset(
    asset_path=Path('/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_Cosmx_Bitmask-scene-000.ome.zarr'),
    label="MERCK_IBD7757A_Cosmx_Bitmask test",
    description="MERCK_IBD7757A_Cosmx_Bitmask tes"
)
bitmask_id = bitmas_up.asset.id #'01922e60-9bf0-8cde-b3c5-3fb6419ae59d'

# IF image
img = Path("/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_Cosmx_IF-scene-000.ome.zarr")
img_up = sdk_client.asset_storage_client.upload_new_asset(
    asset_path=img,
    label="MERCK_IBD7757A_Cosmx_DAPI test",
    description="MERCK_IBD7757A_Cosmx_DAPI test"
)
img_id = img_up.asset.id #01922e72-5420-8500-ae3e-94eaf09b788e

# anndata zarr
adata = Path("/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_Cosmx_raw_anndata.zarr")
adata_up = sdk_client.asset_storage_client.upload_new_asset(
    asset_path=adata,
    label="MERCK_IBD7757A_Cosmx_anndata test",
    description="MERCK_IBD7757A_Cosmx_anndata test"
)
adata_id = adata_up.asset.id #01922e84-224c-835b-bdd6-d54a1581edfe

# anndata zarr 2
adata = Path("/Users/thaotran/GitHub/Project/7_Merck_2024/code/IBD/IBD7757A_Cosmx_merck_anndata.zarr")
adata_up = sdk_client.asset_storage_client.upload_new_asset(
    asset_path=adata,
    label="MERCK_IBD7757A_Cosmx_anndata with cell type test",
    description="MERCK_IBD7757A_Cosmx_anndata with cell type test"
)
adata_id = adata_up.asset.id #01923097-7a08-8a40-afdd-739b4a969221

# tx spot zarr
tx_zarr = Path("/Users/thaotran/GitHub/")
tx_zarr_up = sdk_client.asset_storage_client.upload_new_asset(
    asset_path=tx_zarr,
    label="CosMx Public Half mouse brain tx spot",
    description="CosMx Public Half mouse brain tx spot"
)
tx_zarr_id = tx_zarr_up.asset.id 
