# Output_folder

Written by the scripts; the content is not versioned (see `.gitignore`), only this structure is.

| Folder | Written by | Files (`<tag>` = `run_tag` in the config, e.g. `2025`) |
| --- | --- | --- |
| `01_Preprocessing/` | `Preprocessing_HydraulicModule.m` | `PRE_centerlines_<tag>.csv` (ID, X, Y, z_DTM, W_mask, segment_code, slope), `PRE_segments_endpoints_<tag>.csv`, `PRE_workspace_<tag>.mat`, `figures/PRE_*.png` |
| `02_HydraulicModule/` | `Preprocessing_HydraulicModule.m` | `HM_centerline_riverbed_<tag>.csv`, `HM_h_rect_all_<tag>.mat`, `HM_z_rect_all_<tag>.mat`, `HM_DTM_corrected_<tag>.tif`, `HM_DoD_<tag>.tif`, `HM_workspace_<tag>.mat` (input of the morphodynamic module), `figures/HM_*.png` |
| `03_MorphodynamicModule/` | `MorphodynamicModule.m` | `MM_centerline_riverbed_<tag>.csv`, `MM_h_rect_all_<tag>.mat`, `MM_z_rect_all_<tag>.mat`, `MM_DTM_corrected_<tag>.tif`, `MM_DoD_<tag>.tif`, `MM_workspace_<tag>.mat` (all results, including the node solutions), `figures/MM_*.png` |
| `04_Validation/HydraulicModule/` | `validation_gnss_mean.m` (called by `Preprocessing_HydraulicModule.m`) | `gnss_validation_mean_records.csv`, `gnss_validation_mean_metrics.csv`, `gnss_validation_mean.png` |
| `04_Validation/MorphodynamicModule/` | `validation_gnss_mean.m` (called by `MorphodynamicModule.m`) | same files as above, for the morphodynamic DTM |
| `04_Validation/ModuleComparison/` | `ModuleComparison.m` | `discharge_balance_nodes.csv`, `discharge_balance_segments.csv`, `GNSS_stations_reach_class.csv`, `GNSS_module_comparison_stations.csv`, `GNSS_module_comparison_groups.csv`, `GNSS_module_comparison.png`, `GNSS_station_assignment.png`, `DTM_difference_MM_minus_HM.tif`, `DTM_difference_MM_minus_HM.png`, `DTM_difference_MM_minus_HM_stats.csv`, `GNSS_cross_section_profiles.png`, `GNSS_cross_section_profiles_sections.csv` |
