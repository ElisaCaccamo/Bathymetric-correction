# Output_folder

Written by the scripts; the content is not versioned (see `.gitignore`), only this structure is.

| Folder | Written by | Files (`<tag>` = `run_tag` in the config, e.g. `2025`) |
| --- | --- | --- |
| `01_Preprocessing/` | `Preprocessing_HydraulicModule.m` | `PRE_centerlines_<tag>.csv` (ID, X, Y, z_DTM, W_mask, segment_code, slope), `PRE_segments_endpoints_<tag>.csv`, `PRE_workspace_<tag>.mat`, `figures/PRE_*.png` |
| `02_HydraulicModule/` | `Preprocessing_HydraulicModule.m` | `HM_centerline_riverbed_<tag>.csv`, `HM_h_rect_all_<tag>.mat`, `HM_z_rect_all_<tag>.mat`, `HM_DTM_corrected_<tag>.tif`, `HM_DoD_<tag>.tif`, `HM_workspace_<tag>.mat` (input of the morphodynamic module), `figures/HM_*.png` |
| `03_MorphodynamicModule/` | morphodynamic module (to be adapted) | — |
| `04_Validation/HydraulicModule/` | `validation_gnss_mean.m` | `gnss_validation_mean_records.csv`, `gnss_validation_mean_metrics.csv`, `gnss_validation_mean.png` |

Correspondence with the old names of `0_2025_ORCO_modifiche.mlx`:

| Old | New |
| --- | --- |
| `0_centerline_bat_rect_2025_modifiche_3.csv` | `02_HydraulicModule/HM_centerline_riverbed_2025.csv` |
| `h_rect_all_modifiche_3.mat` | `02_HydraulicModule/HM_h_rect_all_2025.mat` |
| `z_rect_all_modifiche_3.mat` | `02_HydraulicModule/HM_z_rect_all_2025.mat` |
| `0_DTM_bat_corrected_2025_modifiche_5.tif` | `02_HydraulicModule/HM_DTM_corrected_2025.tif` |
| `0_DoD_bat_corrected_2025_modifiche_5.tif` | `02_HydraulicModule/HM_DoD_2025.tif` |
| `hydraulic_workspace_2025.mat` (expected by the morphodynamic module) | `02_HydraulicModule/HM_workspace_2025.mat` |
| `validazione_2025/` | `04_Validation/HydraulicModule/` |
