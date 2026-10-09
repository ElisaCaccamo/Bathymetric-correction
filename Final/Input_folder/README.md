# Input_folder

All the inputs of the pre-processing, hydraulic module, morphodynamic module and validation.
File names and parameters are set in `config_ORCO_2025.m`.

| Folder | Content | File name expected by `config_ORCO_2025.m` |
| --- | --- | --- |
| `DTM_segments/` | DTM of each river segment (filled) | `dtm_000.tif`, `dtm_001.tif`, … |
| `Mask_segments/` | Binary mask of each river segment, same codes as the DTMs | `I_000.tif`, `I_001.tif`, … |
| `DTM/` | Full LiDAR DTM, 50 cm, UTM32N ETRF2000 | `cutDTM_50cm_ORCO_UTM32N_ETRF2000.tif` |
| `WetArea_mask/` | Wet area mask on the same grid as the DTM | `I_wetarea_2025.tif` |
| `GNSS/` | GNSS bed points of 07/03/2025 (`NomePunto, Y, X, Quota ellissoidica, Quota ortometrica`) | `punti_bat.csv` |

`DTM_segments/` and `Mask_segments/` must contain the same number of files, with the same codes:
the i-th DTM is paired with the i-th mask (alphabetical order).
