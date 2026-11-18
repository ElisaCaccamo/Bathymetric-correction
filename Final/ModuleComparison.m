%% MODULE COMPARISON (ORCO)
% Compares the raw LiDAR DTM and the DTMs corrected by the hydraulic (HM) and
% morphodynamic (MM) modules with the GNSS cross-sections, and checks where
% and why the two corrected DTMs differ. Run it after
% Preprocessing_HydraulicModule and MorphodynamicModule (step 3 in README.md).
%
% Terminology: a GNSS "station" is one measured cross-section, i.e. all the
% points named "sNNN_pMMM" with the same "sNNN".
%
% Sections
%   1. Water balance of the flight discharge in both modules (node by node,
%      dead-end branches)                      -> check_discharge_conservation
%   2. Station -> segment -> reach type, with ambiguous assignments and
%      distance from the upstream node         -> station_reach_class
%   3. Statistics by reach type with station-bootstrap intervals, paired
%      test MM vs HM, figure A/B/C, station map -> compare_modules_gnss
%   4. Map of DTM_MM - DTM_HM                   -> dtm_difference_map
%   5. Transverse profiles of selected stations -> gnss_cross_section_profiles
%
% Outputs: Output_folder/04_Validation/ModuleComparison

%% 0. Folders and configuration
root_dir   = fileparts(mfilename('fullpath'));
input_dir  = fullfile(root_dir, 'Input_folder');
output_dir = fullfile(root_dir, 'Output_folder');
config_file = fullfile(input_dir, 'config_ORCO_2025.m');
addpath(root_dir);
run(config_file);

out_hydro  = fullfile(output_dir, '02_HydraulicModule');
out_morpho = fullfile(output_dir, '03_MorphodynamicModule');
out_cmp    = fullfile(output_dir, '04_Validation', 'ModuleComparison');
if ~isfolder(out_cmp), mkdir(out_cmp); end

stations_profiles = {'s006', 's019', 's020'};   % stations of section 5
stations_exclude  = {'s006'};                   % also reported separately in section 3

%% Inputs
ws_hm = load(fullfile(out_hydro, sprintf('HM_workspace_%s.mat', run_tag)), ...
    'all_centerlines', 'sections_ids', 'bif_info', 'conf_info', 'predecessor', 'topo_order', ...
    'network_root', 'Q_map', 'h_rect_all', 'z_rect_all', 'widths_all', 'slope_sect');
ws_mm = load(fullfile(out_morpho, sprintf('MM_workspace_%s.mat', run_tag)), ...
    'Q_map', 'h_rect_morpho', 'z_rect_morpho');
all_centerlines = ws_hm.all_centerlines;

rec_hm = fullfile(output_dir, '04_Validation', 'HydraulicModule', 'gnss_validation_mean_records.csv');
rec_mm = fullfile(output_dir, '04_Validation', 'MorphodynamicModule', 'gnss_validation_mean_records.csv');
T_hm = readtable(rec_hm, 'TextType', 'char');
T_mm = readtable(rec_mm, 'TextType', 'char');

%% 1. Water balance of the flight discharge
[Tnode, Tseg] = check_discharge_conservation(ws_hm.Q_map, ws_mm.Q_map, ws_hm.z_rect_all, ws_mm.z_rect_morpho, ...
    all_centerlines, ws_hm.sections_ids, ws_hm.predecessor, ws_hm.bif_info, ws_hm.conf_info, ...
    ws_hm.network_root, out_cmp);

%% 2. Stations, segments and reach types
S = station_reach_class(T_hm, all_centerlines, ws_hm.sections_ids, ws_hm.predecessor, ...
    ws_hm.bif_info, ws_hm.conf_info, ws_hm.network_root, ws_hm.topo_order);
disp(S);
fprintf('Stations with points on two segments (ambiguous): %s\n', strjoin(S.Station(S.ambiguous)', ', '));
writetable(S, fullfile(out_cmp, 'GNSS_stations_reach_class.csv'));

%% 3. Statistics and comparison figure
[Tgroup, Tst] = compare_modules_gnss(T_hm, T_mm, S, all_centerlines, out_cmp, ...
    struct('exclude', {stations_exclude}, 'B', 2000));

%% 4. Map of the difference between the two corrected DTMs
[dtm_hm, R] = readgeoraster(fullfile(out_hydro,  sprintf('HM_DTM_corrected_%s.tif', run_tag)));
[dtm_mm, ~] = readgeoraster(fullfile(out_morpho, sprintf('MM_DTM_corrected_%s.tif', run_tag)));
[dtm_raw, R_raw] = readgeoraster(fullfile(input_dir, dtm_file));
dtm_hm = double(dtm_hm);  dtm_mm = double(dtm_mm);  dtm_raw = double(dtm_raw);
dtm_raw(dtm_raw == -32767) = NaN;            % NoData of the raw DTM
dtm_hm(dtm_hm <= -9999) = NaN;  dtm_mm(dtm_mm <= -9999) = NaN;
assert(isequal(size(dtm_hm), size(dtm_mm), size(dtm_raw)), 'The three DTMs must be on the same grid.');
[dZ, Tdz] = dtm_difference_map(dtm_hm, dtm_mm, dtm_raw, R, all_centerlines, S, out_cmp, crs_code);

%% 5. Transverse profiles of selected stations
Tsec = gnss_cross_section_profiles(fullfile(input_dir, gnss_file), T_hm, T_mm, stations_profiles, S, ...
    all_centerlines, ws_hm, ws_mm, dtm_raw, R_raw, out_cmp);

disp('Module comparison: all outputs saved in Output_folder/04_Validation/ModuleComparison.');
