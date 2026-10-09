%% CONFIG_ORCO_2025  Parameters of the ORCO 2025 run
% Read by Preprocessing_HydraulicModule.m (and later by the morphodynamic
% module). Everything that was hard-coded or typed with input() in
% 0_2025_ORCO_modifiche.mlx is collected here, so that a run can be repeated
% without interaction (e.g. matlab -batch).

run_tag = '2025';           % used in all output file names

%% Input files (relative to Input_folder)
dtm_segments_pattern  = fullfile('DTM_segments',  'dtm_*.tif');   % DTM of each river segment
mask_segments_pattern = fullfile('Mask_segments', 'I_*.tif');     % binary mask of each river segment
dtm_file      = fullfile('DTM',          'cutDTM_50cm_ORCO_UTM32N_ETRF2000.tif'); % full LiDAR DTM
wetarea_file  = fullfile('WetArea_mask', 'I_wetarea_2025.tif');   % wet area mask (same grid as the DTM)
gnss_file     = fullfile('GNSS',         'punti_bat.csv');        % GNSS points 07/03/2025
crs_code      = 32632;                                            % EPSG of the output GeoTIFFs

%% Pre-processing
plotornot = 0;              % 1 = YES visualization in centerline_from_mask; 0 = NO
Wn = 30;                    % Nominal channel width (average) [m]
dx = 10;                    % Cross-section interval [m]
distance_threshold = 51;    % max distance between end/start of segments to connect them [m]

% River segments to flip UPSTREAM -> DOWNSTREAM (codes as in the mask file names)
segments_to_flip = {'000','006','015','018','019','022','023','024','025', ...
                    '029','030','031','032','033','034','036','041'};

% Bifurcations not found automatically (bridges): one row per bifurcation
% [bifurcating_branch branch1 branch2]
forced_bifurcations = [35 34 36;
                       37 38 39];

% Confluences not found automatically: one row per confluence
% [branch1 branch2 merging_branch]
forced_confluences = zeros(0, 3);

% Targeted smoothing of segments with systematic oscillation of z
segments_to_smooth = struct('tronco', {12, 19, 26, 30, 38, 39, 32}, ...
                            'window', {7,  7,  7,  7,  7,  7,  9});

%% Hydraulic module
g = 9.81;
n_manning = 0.03;           % Manning's coefficient
% LiDAR flight 17/03/2025
Q_upstream = 10.8;          % Discharge SAN BENIGNO ORCO [m^3/s]
level_upstream = 0.695;     % Water level SAN BENIGNO ORCO [m] from most recent (2023) rating curve
slope_min = 1.6e-4;         % minimum slope assigned to a river segment

%% Run options
interactive = false;        % true = ask flips/bifurcations/confluences with input() as in the original
save_figures = true;        % save every figure as PNG in Output_folder
run_diagnostics = true;     % run the diagnostic sections at the end of the hydraulic module
