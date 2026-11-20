%% CONFIG_ORCO_2025  Parameters of the ORCO 2025 run
% Read by Preprocessing_HydraulicModule.m and MorphodynamicModule.m. All the
% parameters of a run are collected here, so that it can be repeated without
% interaction (e.g. matlab -batch).

run_tag = '2025';           % used in all output file names

%% Input files (relative to Input_folder)
dtm_segments_pattern  = fullfile('DTM_segments',  'dtm_*.tif');   % DTM of each river segment
mask_segments_pattern = fullfile('Mask_segments', 'I_*.tif');     % binary mask of each river segment
dtm_file      = fullfile('DTM',          'cutDTM_50cm_ORCO_UTM32N_ETRF2000.tif'); % full LiDAR DTM
wetarea_file  = fullfile('WetArea_mask', 'I_wetarea_2025.tif');   % wet area mask (same grid as the DTM)
gnss_file     = fullfile('GNSS',         'punti_bat.csv');        % GNSS points 07/03/2025
adcp_file     = fullfile('GNSS',         'punti_bat_ADCP.csv');   % ADCP transects 07/03/2025 (Z_BOTTOM ellipsoidal)
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

% Segments without discharge at the flight (codes as in the mask file names).
% They hold still water but carry no appreciable discharge at the flight:
% 15-16 are part of a meander only partly wet, 29 is disconnected from
% segment 22 (connected only at higher discharges). Both modules assign them
% Q = 0 (the whole discharge goes to the other branch of the node, no node
% system), do not compute their depth and leave the raw DTM there.
no_flow_segments = [15 16 29];

% Targeted smoothing of segments with systematic oscillation of z
segments_to_smooth = struct('segment', {12, 19, 26, 30, 38, 39, 32}, ...
                            'window', {7,  7,  7,  7,  7,  7,  9});

%% Validation (field points: GNSS cross-sections by wading + ADCP transects)
% Each ADCP transect is a station (a01, a02, ...) treated like a GNSS cross-section.
% Minimum distance [m] between kept ADCP points along a transect; 0 = all the
% points. STAND-BY: thinning to a spacing similar to the GNSS one (~1.5 m) is a
% question for Carlo (Domande_per_Carlo.md, 6).
adcp_min_spacing = 0;
% GNSS cross-sections that cross more than one segment: split into sub-stations
% (sNNN-tCC), one per segment, each compared with the mean of its own points.
% {station, first point, last point, segment code}; s015 still to be defined.
gnss_station_split = {'s003',  1, 13,  9;
                      's003', 14, 21, 11;
                      's003', 22, 29, 12;
                      's007',  1, 16,  6;
                      's007', 17, 25,  7};

%% Hydraulic module
g = 9.81;
n_manning = 0.03;           % Manning's coefficient
% LiDAR flight 17/03/2025
Q_upstream = 10.8;          % Discharge SAN BENIGNO ORCO [m^3/s]
level_upstream = 0.695;     % Water level SAN BENIGNO ORCO [m] from most recent (2023) rating curve
slope_min = 1.6e-4;         % minimum slope assigned to a river segment

%% Morphodynamic module
% Grain size of the bed surface: sieves [mm] and percent passing of the
% field samples (one row per sample: Pratoregio 28/02/2019, San Benigno 18/03/2019)
grain_sieves_mm = [2 4 8 16 32 64 90 128 256 512 1024 2048];
grain_passing   = [10.80 10.80 10.80 13.20 28.00 56.00 72.40 87.60 100.00 100.00 100.00 100.00;
                    0.25  0.25  0.25  0.25  7.25 47.75 70.25 86.50  99.25 100.00 100.00 100.00];
d_F_rep_mm = 1;             % representative diameter of the fine class (< 2 mm) [mm]
d_ref = 0.001;              % reference diameter of the phi scale [m]
rho = 1000;                 % water density [kg/m^3]
delta = 1.65;               % submerged relative density (rho_s - rho)/rho
teta_r = 0.036;             % Wilcock & Crowe reference Shields number
tau_star_critical = 0.03;   % critical Shields number
bankfull_factor = 1.62;     % bankfull Shields number = bankfull_factor * tau_star_critical
alpha = 4;                  % length of the cells upstream of the node = alpha*W0
r = 0.65;                   % Ikeda coefficient of the transverse gravitational flux
% Slope of the mother at the nodes (fit of the water surface upstream of the node)
L_node = 400;               % window length [m]
L_node_min = 100;           % minimum usable window length [m]
L_node_max = 1600;          % the window is doubled up to this length [m]
n_min_node = 5;             % minimum number of centerline points in the window
S_floor_prctile = 0.05;     % valid slopes are raised to this percentile of the segment slopes
S_fallback_prctile = 0.25;  % slope used when no valid window slope exists
% Bankfull widths on the transect (W_bankfull_from_dtm)
transect_ds = 0.5;          % sampling step [m]
transect_half_len = 100;    % half-length of the transect [m]
spill_band = 3;             % bank tops searched between W_mask/2 and spill_band*W_mask
node_width_source = 'dtm';  % 'dtm' = bankfull widths (method); 'mask' = flight widths (sensitivity test only)
% Node solver
psi_start = [0.10 0.25 0.40 0.60 0.75 0.90];   % starting shares of discharge to branch 1
D_start   = [0.6 1.0 1.4];                     % starting branch depths, multiples of D0
psi_margin = 0.01;          % a physical solution has psi_margin < psi < 1 - psi_margin
tol_res = 1e-8;             % maximum residual norm of an accepted solution

%% Run options
interactive = false;        % true = ask flips/bifurcations/confluences with input()
save_figures = true;        % save every figure as PNG in Output_folder
run_diagnostics = true;     % run the diagnostic sections of the hydraulic and morphodynamic modules
