# Final

Self-contained version of the workflow that reads every input from `Input_folder` and writes every output to `Output_folder`.

| File | Role |
| --- | --- |
| `Preprocessing_HydraulicModule.m` | Pre-processing (network, cross-sections, widths, slopes, topology) + hydraulic module + GNSS validation |
| `MorphodynamicModule.m` | Morphodynamic module (bankfull node model at the bifurcations, fractional sediment transport) + GNSS validation |
| `Input_folder/config_ORCO_2025.m` | All parameters (flips, forced bifurcations, Q, n, sediment and node parameters, ...) |
| `centerline_from_mask_modifiche.m`, `calculate_widths.m`, `skeleton_coords.m`, `crop_to_mask.m` | Functions of the pre-processing |
| `W_bankfull_from_dtm.m`, `G_zi.m` | Functions of the morphodynamic module (bankfull width on the DTM, Wilcock & Crowe transport function) |
| `validation_gnss_mean.m` | GNSS validation of a corrected DTM |

Run from MATLAB with the current folder set to `Final`, in this order:

```matlab
Preprocessing_HydraulicModule
MorphodynamicModule
```

or from a shell: `matlab -batch "Preprocessing_HydraulicModule"`, then `matlab -batch "MorphodynamicModule"`.
See `Input_folder/README.md` for the input files and `Output_folder/README.md` for the outputs.
Theory, hypotheses and equations of the two modules: `docs/TeoriaMetodologica_moduli_idraulico_morfodinamico.md`.
