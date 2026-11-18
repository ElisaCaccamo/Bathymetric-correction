# Final

Self-contained version of the workflow that reads every input from `Input_folder` and writes every output to `Output_folder`.

| File | Role |
| --- | --- |
| `Preprocessing_HydraulicModule.m` | Pre-processing (network, cross-sections, widths, slopes, topology) + hydraulic module + GNSS validation |
| `MorphodynamicModule.m` | Morphodynamic module (bankfull node model at the bifurcations, fractional sediment transport) + GNSS validation |
| `Input_folder/config_ORCO_2025.m` | All parameters (flips, forced bifurcations, Q, n, sediment and node parameters, ...) |
| `centerline_from_mask.m`, `calculate_widths.m`, `skeleton_coords.m`, `crop_to_mask.m` | Functions of the pre-processing |
| `W_bankfull_from_dtm.m`, `G_zi.m` | Functions of the morphodynamic module (bankfull width on the DTM, Wilcock & Crowe transport function) |
| `validation_gnss_mean.m` | GNSS validation of a corrected DTM |
| `ModuleComparison.m` | Comparison of the two modules: water balance of the discharge, GNSS statistics by reach type, map of DTM_MM - DTM_HM, cross-section profiles |
| `check_discharge_conservation.m`, `station_reach_class.m`, `compare_modules_gnss.m`, `dtm_difference_map.m`, `gnss_cross_section_profiles.m` | Functions of `ModuleComparison.m` |

## How to run

The scripts must be run **in this order**, from MATLAB with the current folder set to `Final`:

1. `Preprocessing_HydraulicModule` – builds the river network and runs the hydraulic module. It writes `Output_folder/02_HydraulicModule/HM_workspace_<tag>.mat` and `HM_DTM_corrected_<tag>.tif`.
2. `MorphodynamicModule` – reads those two files, so it can only run after step 1 (it stops with an error if they are missing).
3. `ModuleComparison` (optional) – compares the two modules; it reads the workspaces, the corrected DTMs and the GNSS validation tables of steps 1 and 2, so it runs only after both.

```matlab
Preprocessing_HydraulicModule
MorphodynamicModule
ModuleComparison
```

From a shell: `matlab -batch "Preprocessing_HydraulicModule"`, then `matlab -batch "MorphodynamicModule"`.
Step 1 must be repeated whenever the inputs or the pre-processing/hydraulic parameters in `config_ORCO_2025.m` change; step 2 alone is enough when only the parameters of the morphodynamic module change.

See `Input_folder/README.md` for the input files and `Output_folder/README.md` for the outputs.
Theory, hypotheses and equations of the two modules: `docs/TeoriaMetodologica_moduli_idraulico_morfodinamico.md`.
