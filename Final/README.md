# Final

Self-contained version of the workflow that reads every input from `Input_folder` and writes every output to `Output_folder`.

| File | Role |
| --- | --- |
| `Preprocessing_HydraulicModule.m` | Pre-processing + hydraulic module + GNSS validation (copy of `0_2025_ORCO_modifiche.mlx`) |
| `Input_folder/config_ORCO_2025.m` | All parameters (flips, forced bifurcations, Q, n, ...) |
| `centerline_from_mask_modifiche.m`, `calculate_widths.m`, `skeleton_coords.m`, `crop_to_mask.m`, `validation_gnss_mean.m` | Functions called by the script (copies of the files in the repository root) |

Run from MATLAB with the current folder set to `Final`:

```matlab
Preprocessing_HydraulicModule
```

or from a shell: `matlab -batch "Preprocessing_HydraulicModule"`.
See `Input_folder/README.md` for the input files and `Output_folder/README.md` for the outputs.
