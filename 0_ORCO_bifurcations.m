%% 0_ORCO_bifurcations.m - Morphodynamic module of the bathymetric correction (Orco river)
%
% Reconstructs the bed elevation of the wet area of the Orco river from
% fractional sediment transport at equilibrium and corrects the LiDAR DTM.
% Run 0_2025_ORCO_modifiche.mlx (pre-processing + hydraulic module) first.
% External functions on the path: W_bankfull_from_dtm_modifiche.m, G_zi.m.
% Outputs are saved with the suffix _morpho_BRT.
% Section numbers (§) refer to docs/TeoriaMetodologica_moduli_idraulico_morfodinamico.md.
%
% =========================================================================
% HOW THE MORPHODYNAMIC MODULE WORKS
% =========================================================================
% The river network (centerlines every 10 m, rectangular cross-sections,
% wetted widths W_mask from the wet-area mask of the LiDAR flight, segment
% slopes, topology) comes from the pre-processing. Sections are processed
% from upstream to downstream in the topological order of the segments.
% Each section is a single channel, a bifurcation node or a confluence.
%
% Sediment: two-class mixture, fine F (< 2 mm, d_F = 1 mm) and coarse C
% (d_C = geometric mean of the coarse sieve intervals), surface fractions
% f_F + f_C = 1 (§4.1). Transport per unit width of class x:
% Wilcock & Crowe (2003) with hiding, q_s,x = (tau/rho)^1.5/(g*Delta)*G(zeta_x),
% on the geometric mean diameter d_g of the bed surface (§4.2). The flux of
% class x is Q_s,x = f_x * q_s,x * W.
%
% 1) Single channel (§5.1). Uniform flow (Manning, wide rectangular section)
%    with the discharge from upstream gives the depth D and the bed
%    eta = z_DTM - D. Water discharge Q and sediment fluxes Q_s,F, Q_s,C are
%    conserved along each segment (equilibrium).
%
% 2) Confluence (§5.4). Water and sediment fluxes are summed; the surface
%    composition downstream is the one that, with the local capacity,
%    carries fine and coarse sediment in the ratio of the supply.
%
% 3) Bifurcation node. The split of water and sediment between the two
%    branches is computed at BANKFULL conditions (at the flight discharge
%    the coarse sediment is almost immobile) with the two-cell node model
%    of Bolla Pittaluga, Repetto & Tubino (2003, "BRT"), extended to a
%    two-class mixture as in Ragno et al. (2023):
%    a) Slope of the mother S0: least-squares fit of the water surface over
%       a window L_node = 400 m upstream of the node (local_node_slope).
%    b) Bankfull depth of the mother: Shields criterion
%       D_S = theta_BF*Delta*d_50/S0, theta_BF = 1.62*theta_c, limited by the
%       banks of the transect (spill level = lower of the two bank tops) and
%       never below the critical depth D_c = theta_c*Delta*d_50/S0:
%       D0 = max(min(D_S, z_spill - eta0), D_c) (W_bankfull_from_dtm_modifiche).
%    c) Bankfull width W0 = A0/D0 from the bathymetry-corrected DTM (never
%       below W_mask); branch widths measured in the same way and rescaled
%       so that W0 = W1 + W2. Bankfull discharge Q0 = W0*D0^(5/3)*sqrt(S0)/n.
%       Sediment supply = bankfull capacity of the mother.
%    d) Node system, 6 equations in [psi, D1, D2, S_b, f1F, f2F]
%       (local function morpho_system_brt):
%         - uniform flow in each branch, Q_i = W_i*D_i^(5/3)*sqrt(S_b)/n (2)
%         - cell balance for F and C in two cells of length alpha*W0
%           upstream of the node: lateral sediment flux = flux due to the
%           lateral water flux Q_y minus gravitational flux along the
%           transverse bed slope (Ikeda coefficient r), with
%           eta1 - eta2 = D2 - D1 (common water surface at the node) (2)
%         - sediment continuity for F and C (2)
%       Water continuity Q1 + Q2 = Q0 is built in (Q1 = psi*Q0).
%    e) Choice of the solution. The balanced solution (psi = W1/W0,
%       D1 = D2 = D0, S_b = S0) always exists. Below a critical aspect ratio
%       beta = W0/(2*D0) it is the only solution and it is stable; above it,
%       asymmetric solutions appear and the balanced one is unstable (BRT).
%       Then the asymmetric solution closest to the share of branch 1 in the
%       hydraulic module at flight conditions (psi_hyd) is taken.
%    f) Back to flight conditions: Q1 = psi*Q0_flight, Q2 = (1-psi)*Q0_flight,
%       Q_s1,x = psi_x*Q_s0,x with psi_x = f1x*q_s1x*W1/(f0x*q_s0x*W0);
%       f1F, f2F become the surface composition of the branches. The first
%       section of each branch is computed with uniform flow at its flight
%       discharge, like every other section (eta = z_DTM - D).
%
% 4) DTM correction: the bed elevation of each section is assigned to its
%    cross-section points, interpolated on the DTM grid (scattered
%    interpolant) and written inside the wet-area mask, never above the
%    LiDAR elevation; corrected DTM and DoD are saved as GeoTIFF.
%
% =========================================================================
% HYPOTHESES
% =========================================================================
%  H1  Equilibrium: in every segment the sediment flux of each class equals
%      the transport capacity and is conserved; no erosion/deposition trend.
%  H2  Wide rectangular cross-sections equivalent to the real ones
%      (hydraulic radius = depth; bankfull width = area/depth).
%  H3  Uniform flow (Manning, constant n = 0.03) in every section and in
%      each branch at the node; no backwater, no energy losses at the node.
%  H4  Two grain-size classes (F, C) with the field grain-size curve
%      (Pratoregio and San Benigno, 2019) at the network head; Wilcock &
%      Crowe transport with theta_r = 0.036, Delta = 1.65.
%  H5  The morphology of the nodes is formed at bankfull conditions:
%      theta_BF = 1.62*theta_c, theta_c = 0.03 (source of 1.62 to be
%      documented), limited by the banks and never below theta_c.
%  H6  Common water surface at the node (eta + D equal in the mother and in
%      the two branches) and equal water level at the downstream end of the
%      two branches, which have equal length: hence a common equilibrium
%      slope S_b of the branches (BRT).
%  H7  Two-cell nodal condition: cells of length alpha*W0 (alpha = 4),
%      transverse gravitational flux with Ikeda coefficient r = 0.65 and
%      "grain independence" (theta_x = tau/(rho*g*Delta*d_x)).
%  H8  Branch widths: the DTM bankfull widths fix the ratio W1/W2 and are
%      rescaled so that W0 = W1 + W2.
%  H9  Choice among multiple equilibria: stability (BRT) plus the branch
%      that the hydraulic module indicates as dominant at flight conditions.
%  H10 The bankfull split of water (psi) and sediment (psi_x) also holds at
%      the flight discharge.
%
% =========================================================================
% WHAT IS NEW WITH RESPECT TO EXISTING BIFURCATION MODELS
% =========================================================================
% Reference models: Wang et al. (1995) nodal point relation; Bolla
% Pittaluga, Repetto & Tubino (2003) two-cell model, uniform sediment,
% symmetric branches; Ragno et al. (2021) bifurcation + confluence loop;
% Ragno et al. (2023) two-class mixture, symmetric branches.
%  N1  Inverse use of the node model: the equilibrium configuration is not
%      predicted from a given upstream channel, but used to reconstruct the
%      submerged bed of a real network from LiDAR, wet-area mask and gauging
%      data (bathymetric correction).
%  N2  Branches of different width (W1 ~= W2, W1 + W2 = W0): lateral water
%      flux Q_y and cell balance written with respect to the width-
%      proportional split, so that the balanced solution is psi = W1/W0.
%  N3  Bankfull geometry taken from the DTM at every node (slope window,
%      Shields depth limited by the banks and by the critical depth,
%      equivalent rectangular width), instead of an idealised channel.
%  N4  Choice among multiple equilibria constrained by the observed state
%      (hydraulic module at flight conditions).
%  N5  Coupling with the whole network: water and sediment fluxes and bed
%      surface composition routed through bifurcations and confluences in
%      topological order, then brought back to flight conditions.
% N1-N5 must be checked against the literature before being claimed as
% novel in the paper (in particular N2: extensions of BRT to unequal
% widths may already exist).
%
% =========================================================================
% CHANGES WITH RESPECT TO 2_ORCO_bifurcationsNEW_modifiche_modifiche.mlx
% =========================================================================
%  C1  Node system rewritten in the BRT/Ragno form: uniform flow in the
%      branches added; branch slopes are no longer imposed (common
%      equilibrium slope S_b is an unknown); closure eta0 = (eta1 + eta2)/2
%      removed. Reason: the previous system had no momentum equation in the
%      branches; its asymmetric solutions carried up to 2.5 times the
%      uniform-flow discharge in a branch, and with equal widths it never
%      gave the BRT asymmetric solutions (tested on the 12 nodes, 2025 data).
%  C2  same_slope_node and node_root_rule removed; node_width_source kept
%      ('dtm' default, 'mask' for the sensitivity on beta).
%  C3  Solution choice: stability (balanced if it is the only solution,
%      otherwise asymmetric) and psi_hyd from the hydraulic module.
%  C4  First section of each branch computed with uniform flow at flight
%      conditions (it was set to the bankfull node bed eta_i, which was
%      above the flight water surface at 5-7 nodes).
%  C5  Branch bankfull widths measured with the Shields depth of S0 (the
%      branch slope is an unknown); measured branch slopes kept only as
%      diagnostics.
%  C6  Outputs saved with the suffix _morpho_BRT.
%
% =========================================================================
% TO BE REPORTED IN THE DOCUMENTS AFTER THE RUN IS VERIFIED
% =========================================================================
% docs/TeoriaMetodologica_moduli_idraulico_morfodinamico.md
%  - §4.3: slope of the mother from the window L_node; bankfull depth
%    limited by the spill level and by the critical depth; branch widths
%    with the Shields depth of S0; first branch section by uniform flow.
%  - §5.2-5.3: replace the 8-equation system and its reduction to
%    [q1, eta1, f1F, f2F] with the 6-equation BRT form (C1); remove
%    equation (8); add the stability-based choice of the solution (C3).
%  - §6: update the summary table (bifurcation row).
%  - §7-8: update the open questions and the questions for Carlo with the
%    results of this run (number of balanced / asymmetric nodes, beta,
%    sensitivity to node_width_source, theta_BF).
%  - Change log.
% PaperBAT
%  - Methods 5.3: node system (BRT form), hypotheses H1-H10.
%  - Introduction / contribution: N1-N5 after the literature check; the
%    closure eta0 = (eta1 + eta2)/2 is no longer part of the method.
%  - Results 6.2-6.3 and sensitivity 6.4: validation against GNSS/ADCP by
%    section type; sensitivity to node_width_source (beta), theta_BF/theta_c,
%    r, alpha.
%  - Discussion: equilibrium hypothesis, low banks (critical-depth bound),
%    choice among multiple equilibria.
% Results of the verification run (fill in after the run):
%  - date / LiDAR flight:
%  - nodes balanced / asymmetric / fallback:
%  - GNSS validation (RMSE, bias) hydraulic vs morphodynamic:


%% 1. Workspace and inputs from the hydraulic module
clc;
% Variables produced by the pre-processing and hydraulic module
% (0_2025_ORCO_modifiche.mlx) that this script needs. Everything else in
% the workspace is cleared.
required_vars = {'all_centerlines', 'sections_ids', 'bif_info', 'conf_info', ...
                 'widths_all', 'slope_sect', 'z_rect_all', ...
                 'predecessor', 'topo_order', 'network_root', ...
                 'x_matrix', 'y_matrix', 'z_matrix'};
clearvars('-except', 'required_vars', required_vars{:});

% IMPORT OF VARIABLES FROM THE HYDRAULIC MODULE
%  If 0_2025_ORCO_modifiche.mlx has just been run in the same session, the
%  variables are already in the workspace. Otherwise they are loaded from
%  hydro_file, which is created at the end of the hydraulic script with:
%     save(hydro_file, 'all_centerlines', 'sections_ids', 'bif_info', 'conf_info', ...
%          'widths_all', 'slope_sect', 'z_rect_all', 'predecessor', 'topo_order', ...
%          'network_root', 'x_matrix', 'y_matrix', 'z_matrix');
%  All variables must come from the SAME run of the hydraulic module (same
%  order and flips of the river segments), otherwise z_rect_all and
%  all_centerlines are not aligned.
%
%  all_centerlines columns: [ID, X, Y, z_DTM (water surface), W mask, segment code]
%  slope_sect: segment slope from the pre-processing (linear fit of z_DTM)
%  z_rect_all: riverbed elevation from the hydraulic module
hydro_file = '\\nas-camporeale\home\ELISA\BATHIMETRY\0_MODIFICHE\hydraulic_workspace_2025.mat';

missing_vars = {};
for k = 1:numel(required_vars)
    if ~exist(required_vars{k}, 'var')
        missing_vars{end+1} = required_vars{k}; %#ok<SAGROW>
    end
end
if ~isempty(missing_vars)
    if ~isfile(hydro_file)
        error(['Missing variables from the hydraulic module (%s) and file %s not found. ' ...
               'Run 0_2025_ORCO_modifiche.mlx and save its variables in hydro_file.'], ...
               strjoin(missing_vars, ', '), hydro_file);
    end
    load(hydro_file, missing_vars{:});
end

% Consistency checks on the imported variables
N_chk = size(all_centerlines, 1);
assert(isequal(all_centerlines(:,1), (1:N_chk)'), ...
    'Section IDs (all_centerlines(:,1)) must coincide with the row indices 1..N.');
assert(numel(widths_all) == N_chk && numel(slope_sect) == N_chk && numel(z_rect_all) == N_chk, ...
    'widths_all, slope_sect and z_rect_all do not match all_centerlines: different hydraulic runs?');
assert(size(x_matrix,1) == N_chk && size(y_matrix,1) == N_chk && size(z_matrix,1) == N_chk, ...
    'x_matrix, y_matrix and z_matrix must have one row per section.');
assert(numel(predecessor) == size(sections_ids,1) && ...
       isequal(sort(topo_order(:)), (1:size(sections_ids,1))'), ...
    'predecessor/topo_order are not consistent with sections_ids.');
widths_all = widths_all(:);
slope_sect = slope_sect(:);
z_rect_all = z_rect_all(:);


%% 2. Grain size distribution of the surface mixture (§4.1)
% Global mean grain size distribution of the surface mixture (methodology §4.1)
%  Field samples 02/28/2019 (Pratoregio) and 03/18/2019 (San Benigno meander)
diametri = [2 4	8 16 32	64 90 128 256 512 1024 2048]; % sieve sizes [mm]
passanteA1 = [10.80	10.80	10.80	13.20	28.00	56.00	72.40	87.60	100.00	100.00	100.00	100.00];
passanteA2 = [0.25	0.25	0.25	0.25	7.25	47.75	70.25	86.50	99.25	100.00	100.00	100.00];

passante_media = mean([passanteA1; passanteA2]);   % mean percent passing

% Fine class (<= 2 mm) and coarse class (> 2 mm)
indice_fine = diametri <= 2;
indice_coarse = diametri > 2;

% Fraction of the mixture in each sieve interval
classi_frazioni = diff([0, passante_media]) / 100;

% Fractions of the two classes (upstream mother channel: f_0F, f_0C)
f0f = sum(classi_frazioni(indice_fine));
f0c = sum(classi_frazioni(indice_coarse));

% Representative class diameters in the phi scale (methodology §4.1):
% each sieve interval is represented by its phi midpoint, i.e. by the
% geometric mean diameter sqrt(d_{k-1}*d_k), and d_C = 2^(-sum(f_k*phi_k)/f_C).
% The fine class (< 2 mm) has no lower sieve: d_F = 1 mm (phi = 0).
d_F_rep_mm = 1;
d_mid = sqrt(diametri(1:end-1) .* diametri(2:end));   % [mm] intervals of the coarse class
phi_mid = -log2(d_mid);
phi_C = sum(classi_frazioni(indice_coarse) .* phi_mid) / f0c;
d_c = 2^(-phi_C) / 1000;      % [m]
d_f = d_F_rep_mm / 1000;      % [m]

% d_50 of the mean curve (used only for the bankfull depth, methodology §4.3)
[passante_ordinato, idx_unique] = unique(passante_media, 'stable');
diametri_ordinati = diametri(idx_unique);
d_50 = interp1(passante_ordinato, diametri_ordinati, 50) / 1000;   % [m]

% Check
somma_frazioni = f0f + f0c;
fprintf('Fine sediment fraction (<2 mm): %.2f %%\n', 100*f0f);
fprintf('Coarse sediment fraction (>=2 mm): %.2f %%\n', 100*f0c);
fprintf('Sum: %.2f %%\n', 100*somma_frazioni);
fprintf('d_F = %.1f mm, d_C = %.1f mm, d_50 = %.1f mm\n', 1000*d_f, 1000*d_c, 1000*d_50);

figure;
semilogx(diametri, passante_media, '-o', 'LineWidth', 2);
grid on;
xlabel('Diameter (mm)', 'FontSize', 12);
ylabel('Passing (%)', 'FontSize', 12);
title('Mean grainsize distribution', 'FontSize', 14);
xlim([1 3000]);
ylim([0 100]);
hold on;
% Show d_50 in the plot
plot(d_50*1000, 50, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
text(d_50*1000*1.1, 50, sprintf('d_{50} = %.1f mm', d_50*1000), 'FontSize', 12);


%% 3. Bankfull width from the bathymetry-corrected DTM (§4.3)
% Folders with the input rasters (wet-area mask, LiDAR DTM) and with the
% functions G_zi.m and W_bankfull_from_dtm_modifiche.m
addpath('\\nas-camporeale\home\ELISA\BATHIMETRY\2025_bat');
addpath('\\nas-camporeale\home\ELISA\BATHIMETRY\0_MODIFICHE');
% The morphodynamic approach works at bankfull conditions (methodology §4.3):
% at the discharge of the LiDAR flight the coarse sediment is almost
% immobile. A bankfull width corresponding to the bankfull depth is needed.

% W_bankfull_from_dtm_modifiche (external function) returns the theoretical
% bankfull width of a cross-section:
% • Bankfull level z_target = z_rect_all(i) + D, measured from the bed of the
% hydraulic module (D = Shields depth with the slope S0 of the mother, for
% the mother section and for the first section of each branch).
% • Transect orthogonal to the centerline, with the direction computed only
% from points of the same river segment, from -100 to +100 m every 0.5 m; the
% corrected DTM is sampled along it.
% • Banks: first crossings of the profile with z_target on each side, outside
% the wetted channel (|s| >= W_mask/2). If only one bank is found, the other
% one is mirrored.
% • Area A between the banks and equivalent rectangular width W_bf = A /
% (z_target - z_rect_all), consistent with the rectangular schematization (A =
% D * W).
% • Bankfull level limited by the banks (§7.3): bank tops = highest terrain
% point on each side between W_mask/2 and 3*W_mask; spill level z_spill = the
% lower of the two. If the Shields level is higher, D_bf = z_spill -
% z_rect_all (1 cm below z_spill, so both banks are crossed) and the function
% returns this D_bf, which is then used for W, Q_bf, the bankfull transport
% and the node system. This removes most clamps and mirrored banks. Lower
% bound: the returned depth is never below the critical depth
% theta_c*delta*d_50/S (low banks, e.g. a bar close to the channel); in that
% case the width stays the one measured at the spill level and the section is
% flagged (level code 2).
% • If W_bf < W_mask (wetted width of the LiDAR flight), W_bf = W_mask
% (opts.W_old_rule = 'max'; 'nan' discards the section, for debugging).

%% --- Load the bathymetry-corrected DTM (output of the hydraulic module) ---
% Used ONLY to measure bankfull widths (W_bankfull_from_dtm_modifiche): inside the wet
% area its bed is z_rect_all, consistent with the bankfull level
% z_rect_all + D and with the area integrated between the banks; outside
% the channel it coincides with the LiDAR DTM.
dtm_file = '\\nas-camporeale\home\ELISA\BATHIMETRY\0_MODIFICHE\0_DTM_bat_corrected_2025_modifiche_5.tif';
[Z, R] = readgeoraster(dtm_file);   % Z = matrix, R = referencing object
Z = double(Z);
% NoData: read from the file (raw DTM uses -32767, corrected DTM uses NaN);
% clearly non-physical values are discarded anyway.
info_dtm = georasterinfo(dtm_file);
nodata_val = info_dtm.MissingDataIndicator;
if ~isempty(nodata_val)
    Z(Z == nodata_val) = NaN;
end
Z(Z <= -9999) = NaN;

% World coordinates of the columns (Xgrid) and rows (Ygrid) of Z.
% intrinsicToWorld works for both cell rasters (MapCellsReference) and
% posting rasters (MapPostingsReference), and for rows stored from north
% or from south, so Xgrid/Ygrid always match the order of Z for interp2.
[nRows, nCols] = size(Z);
[Xgrid, ~] = R.intrinsicToWorld(1:nCols, ones(1, nCols));
[~, Ygrid] = R.intrinsicToWorld(ones(1, nRows), 1:nRows);


%% 4. Morphodynamic module: single channel, bifurcation node and confluence (§5)
% • Single channel (§5.1): uniform-flow depth D, bed eta = z_DTM - D; sediment
% fluxes Q_s,x and surface composition constant along each segment.
% • Bifurcation (see the header, point 3): slope of the mother S0 from a
% window upstream of the node; bankfull depth D0 (Shields, limited by the
% banks and by the critical depth), width W0 from the DTM, Q_bf from Manning;
% branch widths from the DTM rescaled so that W1 + W2 = W0; node system in
% the BRT form, 6 equations in [psi, D1, D2, S_b, f1F, f2F] (uniform flow in
% the branches, cell balance and sediment continuity for F and C), solved
% from several starting points; balanced solution if it is the only one,
% otherwise the asymmetric solution closest to the hydraulic share psi_hyd;
% back to flight conditions with psi (water) and psi_x (sediment per class).
% • Confluence (§5.4): discharges and sediment fluxes summed; downstream
% surface composition in equilibrium with the supply.
%% ======================= PARAMETERS =======================
rho = 1000;       % water density [kg/m^3]
g = 9.81;         % [m/s^2]
n_manning = 0.03; % Manning coefficient (k_s = 1/n, same as the hydraulic module)
teta_r = 0.036;   % Wilcock & Crowe reference Shields number (gravel mixture, methodology §4.2)

tau_star_critical = 0.03;
tau_star = 1.62 * tau_star_critical; % bankfull Shields number, single-thread gravel-bed rivers (§4.3)
delta = 1.65;     % submerged relative density (rho_s - rho)/rho

alpha = 4;       % length of the upstream cells = alpha*W0 (Ragno et al. 2021, 2023)
r = 0.65;        % Ikeda coefficient [0.3-1] (Ikeda 1982; Talmon et al. 1995)
d_ref = 0.001;   % [m] reference grain size of the phi scale (Ragno et al. 2023)

% phi scale of the two classes (§4.1). The geometric mean diameter d_g
% depends on the surface composition f_F and is computed inside
% local_qs_fine_coarse / local_dg.
fic = -log2(d_c / d_ref);
fif = -log2(d_f / d_ref);
fim = f0c * fic + f0f * fif;
sigma = sqrt(f0c * (fic - fim)^2 + f0f * (fif - fim)^2);
% ks_c, ks_f: roughness HEIGHTS [m] (Ragno et al. 2023), reported for
% reference only; uniform flow uses k_s = 1/n_manning.
ks_c = 2 * d_c * 2^sigma;
ks_f = 2 * d_f * 2^sigma;

% Discharge at the network head (San Benigno gauge) during the LiDAR flight
% LiDAR 28/03/2023: Q = 1.68 m^3/s, h = 0.615 m
% LiDAR 16/05/2024: Q = 11.3 m^3/s, h = 0.70 m
% LiDAR 17/03/2025: Q = 10.8 m^3/s, h = 0.695 m
Q_upstream = 10.8;

N = size(all_centerlines, 1);

section_to_index = @(id) id;   % section ID = row index (checked in the import block)

%% ======================= PREALLOCATION =======================
% --- Bankfull quantities (bifurcation nodes and their branches) ---
Q_bf_all = NaN(N,1);
h_bf_all = NaN(N,1);           % mother: D0 (Shields, limited by the banks); branches: D1, D2 from the node system
W_bf_all = NaN(N,1);           % mother: from the DTM; branches: W0 split with the DTM ratio W1^bf/W2^bf
W_bf_dtm_all = NaN(N,1);       % branches: bankfull width from the DTM with the Shields depth of S0 (before rescaling)
A_bf_all = NaN(N,1);
left_pt_all = NaN(N,2); right_pt_all = NaN(N,2);
qsf_bf_all = NaN(N,1); qsc_bf_all = NaN(N,1);   % bankfull transport capacity per unit width
f1f_all = NaN(N,1); f1c_all = NaN(N,1);         % surface composition of branch 1 (node solution)
f2f_all = NaN(N,1); f2c_all = NaN(N,1);         % surface composition of branch 2 (node solution)
q1_all = NaN(N,1); eta1_all = NaN(N,1);         % bankfull solution: q per unit width [m^2/s], bed elevation eta0 + D0 - D1 [m] (diagnostics)
q2_all = NaN(N,1); eta2_all = NaN(N,1);
h1_bf_all = NaN(N,1); h2_bf_all = NaN(N,1);     % bankfull depths D1, D2 of the branches
bif_exitflag_all = NaN(N,1);                    % fsolve exit flag at each mother section (diagnostics)
bif_resnorm_all = NaN(N,1);                     % norm of the node residual at the solution (diagnostics)
S_node_all = NaN(N,1);                          % mother row -> S0 used in the node; first row of each branch -> slope MEASURED downstream (diagnostics only)
Sb_all = NaN(N,1);                              % equilibrium slope of the branches from the node system (mother row)
psi_hyd_all = NaN(N,1);                         % share of branch 1 in the hydraulic module at flight conditions (mother row)
S_node_src_all = NaN(N,1);                      % origin of S_node_all (codes in local_node_slope)
D_shields_all = NaN(N,1);                       % Shields bankfull depth tau_star*delta*d_50/S (before the bank limit)
level_code_all = NaN(N,1);                      % bankfull level: 0 = Shields, 1 = spill level, 2 = spill + critical depth (low banks)
node_solutions = cell(N,1);                     % all converged solutions of each node (diagnostics, see cell 5)
beta_all = NaN(N,1);                            % aspect ratio of the mother at bankfull, W0/(2*D0)
theta_bf_all = NaN(N,1);                        % bankfull Shields number of the mother D0*S0/(delta*d_50)
psi_all = NaN(N,1);                             % bankfull share of the water discharge to branch 1 (selected solution)
bif_type_all = NaN(N,1);                        % 0 = fallback, 1 = balanced, 2 = asymmetric solution
bif_nsol_all = NaN(N,1);                        % number of distinct physical solutions found
% --- Quantities at the conditions of the LiDAR flight ---
W_norm_all = NaN(N,1);                          % branch widths brought back to flight conditions
q_norm_all = NaN(N,1);                          % branch discharge [m^3/s] at flight conditions
qsf_cap_all = NaN(N,1); qsc_cap_all = NaN(N,1); % Wilcock-Crowe capacity per unit width
qs0f_all = NaN(N,1); qs0c_all = NaN(N,1);       % actual transport per unit width (Q_s,x / W)
qs_tot_f_all = NaN(N,1);                        % total sediment flux Q_s,x [m^3/s]
qs_tot_c_all = NaN(N,1);
total_QS_all_f = NaN(N,1);
total_QS_all_c = NaN(N,1);
f0f_all = NaN(N,1);                             % surface composition of each section
f0c_all = NaN(N,1);

% Transect parameters for W_bankfull_from_dtm
opts = struct();
opts.ds = 0.5;                 % sampling step along the transect [m]
opts.max_half_len = 100;       % transect half-length [m]
opts.interp_method = 'linear';
opts.h_bf_all = h_bf_all;
opts.W_old_rule = 'max';       % W_bf < W mask -> W_bf = W mask ('nan' discards the section, for debugging)
opts.verbose = false;          % true to print diagnostics for each transect
opts.limit_to_spill = true;    % bankfull level limited to the spill level of the transect (§7.3)
opts.spill_band = 3;           % bank tops searched between W_mask/2 and spill_band*W_mask

h_rect_morpho = NaN(N,1);      % depth at flight conditions
z_rect_morpho = NaN(N,1);      % bed elevation of the morphodynamic module

sectionWet(N) = struct('h', NaN, 'quot', NaN, 'Q', NaN);
hydro_data(N) = struct('W_bf',NaN,'A_bf',NaN,'left_pt',[NaN NaN],'right_pt',[NaN NaN],...
                       'h_bf',NaN,'Q_bf',NaN,'W_norm',NaN,'q_norm',NaN,'h_morph',NaN,'z_morph',NaN);

ks_cd50 = []; % optional Strickler coefficient [m^(1/3)/s]: [] = 1/n, scalar or [N x 1] vector
ks_i_default = 1 / n_manning;

sediment_data = struct( ...
    'qs_c', num2cell(NaN(N,1)), ...
    'qs_f', num2cell(NaN(N,1)), ...
    'fc',   num2cell(NaN(N,1)), ...
    'ff',   num2cell(NaN(N,1)), ...
    'Qs_tot_f',   num2cell(NaN(N,1)), ...
    'Qs_tot_c',   num2cell(NaN(N,1)));

%% ======================= NETWORK TOPOLOGY =======================
bif_ids  = cell2mat_safe(bif_info);      % mother sections (last point of the mother segment)
conf_ids = cell2mat_safe(conf_info);     % confluence sections (first point of the downstream segment)
bif_child_ids = [];                      % first point of each bifurcation branch
for b = 1:size(bif_info,1)
    bif_child_ids = [bif_child_ids, bif_info{b,2}(:)']; %#ok<AGROW>
end
seg_of_row = zeros(N,1);                 % segment index of each section
for s = 1:size(sections_ids,1)
    seg_of_row(sections_ids(s,1):sections_ids(s,2)) = s;
end
% Downstream segments of each segment (inverse of predecessor), used to
% extend the slope windows of the branches across segment boundaries.
successor = cell(size(sections_ids,1), 1);
for s = 1:size(sections_ids,1)
    for u = reshape(predecessor{s}, 1, [])
        successor{u}(end+1) = s;
    end
end
is_segment_start = false(N,1);
is_segment_start(sections_ids(:,1)) = true;

% Water discharge: assigned at the network head (as in the hydraulic
% module) and then always taken from upstream (previous section of the
% same segment, end of the upstream segment, split at a bifurcation, sum
% at a confluence).
Q_map = NaN(N,1);
Q_map(sections_ids(network_root,1)) = Q_upstream;

% Surface composition and sediment fluxes per segment (§5.1, §5.4):
% f_F = f_0F at the network head, inherited along simple continuations,
% f_1F / f_2F in the bifurcation branches, recomputed downstream of each
% confluence. Q_s,x is constant along each segment (equilibrium).
n_seg = size(sections_ids,1);
f_seg_f = NaN(n_seg,1);
f_seg_f(network_root) = f0f;
Qs_seg = NaN(n_seg,2);          % [Q_s,F  Q_s,C] per segment [m^3/s]

% Processing order: segments in topological order (topo_order of the
% hydraulic module), points of each segment from upstream to downstream.
proc_order = cell2mat(arrayfun(@(s) (sections_ids(s,1):sections_ids(s,2))', ...
                      topo_order(:), 'UniformOutput', false));

% ----- Slope of the mother at the bifurcation nodes -----
% The segment slopes slope_sect come from a least-squares fit over the WHOLE
% segment and are floored at 1.6e-4 in the pre-processing: segments with a
% flat or reversed fit get S = 1.6e-4, i.e. a Shields bankfull depth of
% ~31 m. At each node the slope of the mother S0 is therefore recomputed on
% the water surface (all_centerlines(:,4)) over a window of length L_node
% upstream of the node, crossing segment boundaries when the segment is
% shorter (local_node_slope). The slopes of the branches are NOT data of
% the node: they are unknowns of the node system (equilibrium slope S_b, see
% below). The slopes measured downstream of the node are computed with the
% same window only for comparison (diagnostics).
L_node = 400;          % [m] window length (300-500 m, i.e. ~10 times the mask width of the mothers)
L_node_max = 1600;     % [m] the window is doubled up to this length when the fit is not usable
L_node_min = 100;      % [m] minimum usable window length (nodes close to the network head/outlet)
n_min_node = 5;        % minimum number of centerline points in the window
% Reference slopes from the data (instead of the 1.6e-4 floor of the
% pre-processing), computed on the segment slopes that are not at that
% floor. slope_sect has one value per section, so the percentiles are
% weighted by segment length.
%  S_floor_node   : 5th percentile. A valid window slope is kept as measured
%                   and only raised to S_floor_node (values close to zero).
%  S_fallback_node: 25th percentile. Used only when no valid window slope
%                   can be computed.
% (Using the 25th percentile as a floor for every node, as in the previous
% run, overwrote 9 of the 12 measured mother slopes.)
S_floor_pre = 1.6e-4;
S_valid = sort(slope_sect(isfinite(slope_sect) & slope_sect > 1.001 * S_floor_pre));
if isempty(S_valid)
    S_floor_node = S_floor_pre;
    S_fallback_node = S_floor_pre;
else
    S_floor_node    = S_valid(max(1, ceil(0.05 * numel(S_valid))));
    S_fallback_node = S_valid(max(1, ceil(0.25 * numel(S_valid))));
end
fprintf('Node slopes: window L = %.0f m (min %.0f m, max %.0f m), floor = %.2e, fallback = %.2e\n', ...
    L_node, L_node_min, L_node_max, S_floor_node, S_fallback_node);

% Widths used in the node:
%   'dtm' : bankfull widths measured on the bathymetry-corrected DTM
%           (W_bankfull_from_dtm_modifiche), branch widths rescaled on W0;
%   'mask': widths of the wet-area mask of the LiDAR flight (sensitivity test:
%           the aspect ratio beta = W0/(2*D0), which controls whether
%           asymmetric solutions exist, is then much lower).
node_width_source = 'dtm';

% Node solver: fsolve is started from the balanced configuration and from
% asymmetric configurations (share of branch 1 psi_start, branch depths
% D_start*D0), so that every solution of the node system is found.
psi_start = [0.10 0.25 0.40 0.60 0.75 0.90];
D_start   = [0.6 1.0 1.4];
psi_margin = 0.01;     % a physical solution has psi_margin < psi < 1 - psi_margin
tol_res = 1e-8;        % maximum residual norm of an accepted solution

options = optimoptions('fsolve','Display','off','TolFun',1e-10,'TolX',1e-10,'MaxIter',400);

%% ======================= LOOP OVER THE SECTIONS =======================
for kk = 1:numel(proc_order)
    i = proc_order(kk);
    id_i = all_centerlines(i,1);
    s_i = seg_of_row(i);
    is_bif   = ismember(id_i, bif_ids);
    is_conf  = ismember(id_i, conf_ids);
    is_child = ismember(id_i, bif_child_ids);

    % ==================== A) DISCHARGE AND COMPOSITION FROM UPSTREAM ====================
    if ~is_conf && isnan(Q_map(i))
        if is_segment_start(i)
            up = predecessor{s_i};
            if is_child
                warning('Section %d: bifurcation branch without discharge (mother not computed) -> skip', id_i);
                continue
            elseif isempty(up)
                warning(['Segment %d: network head other than network_root, upstream discharge ' ...
                         'not defined -> skip'], s_i);
                continue
            elseif numel(up) == 1
                Q_map(i) = Q_map(sections_ids(up,2));
            else
                error('Segment %d: several upstream segments but section %d is not in conf_info.', s_i, id_i);
            end
        else
            Q_map(i) = Q_map(i-1);   % next section of the same segment
        end
    end
    if is_segment_start(i) && ~is_conf && ~is_child && isnan(f_seg_f(s_i))
        up = predecessor{s_i};
        if numel(up) == 1   % simple continuation: composition and sediment fluxes unchanged
            f_seg_f(s_i) = f_seg_f(up);
            Qs_seg(s_i,:) = Qs_seg(up,:);
        end
    end

    % ==================== B) HYDRAULICS OF THE SECTION ====================
    W_i = widths_all(i);
    if is_conf
        % ---------------- CONFLUENCE (§5.4) ----------------
        conf_idx = find(conf_ids == id_i);
        upstream_sections = conf_info{conf_idx, 2};
        idx1 = section_to_index(upstream_sections(1));
        idx2 = section_to_index(upstream_sections(2));
        Q1 = Q_map(idx1);
        Q2 = Q_map(idx2);
        if isnan(Q1) && isnan(Q2)
            warning('Confluence %d: discharge not defined in both branches -> skip', id_i);
            continue
        elseif isnan(Q1) || isnan(Q2)
            warning('Confluence %d: discharge missing in one branch, using the other one only', id_i);
        end
        Q_map(i) = sum([Q1, Q2], 'omitnan');
        h_i = local_normal_h(Q_map(i), i, widths_all, slope_sect, ...
            all_centerlines, h_rect_morpho, ks_cd50, ks_i_default);
        if ~isfinite(h_i) || h_i <= 0
            warning('Confluence %d: non-physical depth (%.3g) -> skip', id_i, h_i);
            continue
        end

        % Sediment fluxes: sum of the inputs (Q_s,x constant along the branches)
        Qs_up = Qs_seg(seg_of_row([idx1 idx2]), :);
        Qs_seg(s_i,:) = sum(Qs_up, 1, 'omitnan');

        % Downstream surface composition: the f_F that, with the local
        % Wilcock-Crowe capacity (which depends on f through d_g), carries
        % fine and coarse sediment in the ratio of the supply:
        %   f*qsF(f) / ((1-f)*qsC(f)) = Q_s,F / Q_s,C
        f_up = f_seg_f(seg_of_row([idx1 idx2]));
        f_seg_f(s_i) = local_equilibrium_f(Qs_seg(s_i,1), Qs_seg(s_i,2), rho, g, ...
            slope_sect(i), h_i, d_f, d_c, d_ref, teta_r, delta, mean(f_up, 'omitnan'));
    else
        % ---------------- SINGLE CHANNEL / MOTHER SECTION / FIRST SECTION OF A BRANCH ----------------
        % Uniform flow with the discharge from upstream. For the first
        % section of a branch the discharge is psi*Q0 (or (1-psi)*Q0) from
        % the node system: the morphodynamic module acts on the bed through
        % the split of water and sediment, while the bed at flight conditions
        % always comes from uniform flow, eta = z_DTM - D, as in every other
        % section.
        h_i = local_normal_h(Q_map(i), i, widths_all, slope_sect, ...
            all_centerlines, h_rect_morpho, ks_cd50, ks_i_default);
        if ~isfinite(h_i) || h_i <= 0
            warning('Section %d: non-physical depth (%.3g) -> skip', i, h_i);
            continue;
        end
    end

    h_rect_morpho(i) = h_i;
    z_rect_morpho(i) = all_centerlines(i,4) - h_i;     % eta = z_DTM - D
    sectionWet(i).h = h_i;
    sectionWet(i).quot = z_rect_morpho(i);
    sectionWet(i).Q = Q_map(i);
    hydro_data(i).W_norm = W_i;
    hydro_data(i).q_norm = Q_map(i);
    hydro_data(i).h_morph = h_i;
    hydro_data(i).z_morph = z_rect_morpho(i);

    % ==================== C) SEDIMENT TRANSPORT OF THE SECTION ====================
    f_i_f = f_seg_f(s_i);
    f_i_c = 1 - f_i_f;
    [qsf_cap, qsc_cap] = local_qs_fine_coarse(rho, g, slope_sect(i), h_i, ...
        d_f, d_c, d_ref, f_i_f, teta_r, delta);
    qsf_cap_all(i) = qsf_cap;
    qsc_cap_all(i) = qsc_cap;
    if any(isnan(Qs_seg(s_i,:)))
        % first section of a segment with no known supply (network head):
        % Q_s,x = f_x * q_s,x * W
        Qs_seg(s_i,:) = [f_i_f * qsf_cap, f_i_c * qsc_cap] * W_i;
    end
    qs_tot_f_all(i) = Qs_seg(s_i,1);
    qs_tot_c_all(i) = Qs_seg(s_i,2);
    total_QS_all_f(i) = qs_tot_f_all(i);
    total_QS_all_c(i) = qs_tot_c_all(i);
    qs0f_all(i) = qs_tot_f_all(i) / W_i;   % changes with the width, the total flux does not
    qs0c_all(i) = qs_tot_c_all(i) / W_i;
    f0f_all(i) = f_i_f;
    f0c_all(i) = f_i_c;
    sediment_data(i).qs_f = qs0f_all(i);
    sediment_data(i).qs_c = qs0c_all(i);
    sediment_data(i).ff = f_i_f;
    sediment_data(i).fc = f_i_c;
    sediment_data(i).Qs_tot_f = qs_tot_f_all(i);
    sediment_data(i).Qs_tot_c = qs_tot_c_all(i);

    if ~is_bif
        continue
    end

    % ==================== D) BIFURCATION NODE (§4.3, §5.2-5.3) ====================
    section_mother = id_i;
    row_bif = find(bif_ids == section_mother);
    figlie = bif_info{row_bif, 2};
    id1 = section_to_index(figlie(1));
    id2 = section_to_index(figlie(2));
    s1 = seg_of_row(id1);
    s2 = seg_of_row(id2);

    % --- Slopes: S0 of the mother (window upstream of the node); the slopes
    % measured downstream of the node are stored for comparison only ---
    net = struct('all_centerlines', all_centerlines, 'sections_ids', sections_ids, ...
                 'seg_of_row', seg_of_row, 'predecessor', {predecessor}, ...
                 'successor', {successor}, 'widths_all', widths_all, ...
                 'Q_map', Q_map, 'slope_sect', slope_sect);
    [S0, src0] = local_node_slope(i,   -1, L_node, L_node_min, L_node_max, n_min_node, S_floor_node, S_fallback_node, net);
    [S1m, src1] = local_node_slope(id1, +1, L_node, L_node_min, L_node_max, n_min_node, S_floor_node, S_fallback_node, net);
    [S2m, src2] = local_node_slope(id2, +1, L_node, L_node_min, L_node_max, n_min_node, S_floor_node, S_fallback_node, net);
    S_node_all([i id1 id2]) = [S0 S1m S2m];
    S_node_src_all([i id1 id2]) = [src0 src1 src2];

    % --- Bankfull conditions of the mother (§4.3) ---
    % Shields depth D_S = tau_star*delta*d_50/S0, limited by the banks inside
    % W_bankfull_from_dtm_modifiche (D0 = min(D_S, z_spill - eta0)) and never
    % below the critical depth D_c = theta_c*delta*d_50/S0 (sediment at least
    % at incipient motion). The SAME D0 is used for W0, Q_bf, the bankfull
    % transport of the mother and the node system.
    D0_S = (tau_star * delta * d_50) / S0;
    opts.h_bf_all = h_bf_all;
    opts.h_min = (tau_star_critical * delta * d_50) / S0;
    [W0, left_pt, right_pt, A_bf, D0, code0] = W_bankfull_from_dtm_modifiche( ...
        i, z_rect_all, D0_S, Z, Xgrid, Ygrid, all_centerlines, opts);   % W0 = A0/D0 (>= W mask)
    if ~(isfinite(D0) && D0 > 0)
        D0 = D0_S;
    end
    if strcmp(node_width_source, 'mask') && isfinite(widths_all(i)) && widths_all(i) > 0
        W0 = widths_all(i);          % sensitivity test: width of the flight mask
    end
    D_shields_all(i) = D0_S;
    level_code_all(i) = code0;
    beta_all(i) = W0 / (2 * D0);
    theta_bf_all(i) = D0 * S0 / (delta * d_50);
    Q_bf = (1 / n_manning) * W0 * D0^(5/3) * sqrt(S0);   % uniform flow of the mother at bankfull
    if ~isfinite(Q_bf) || Q_bf <= 0
        warning('Section %d: invalid Q_bf (%.3f)', i, Q_bf);
        Q_bf = NaN;
    end
    W_bf_all(i) = W0; A_bf_all(i) = A_bf; h_bf_all(i) = D0; Q_bf_all(i) = Q_bf;
    left_pt_all(i,:) = left_pt; right_pt_all(i,:) = right_pt;
    hydro_data(i).W_bf = W0; hydro_data(i).A_bf = A_bf; hydro_data(i).h_bf = D0;
    hydro_data(i).Q_bf = Q_bf; hydro_data(i).left_pt = left_pt; hydro_data(i).right_pt = right_pt;

    % --- Bankfull widths of the branches (§4.3) ---
    % The bankfull width of each branch is measured on the DTM at the level
    % of the Shields depth computed with S0 (the slope of the branches is an
    % unknown of the node, so the mother slope is the only data available),
    % with the same bank limit and critical-depth bound as the mother. The
    % widths are rescaled on W0 to close the node, W0 = W1 + W2.
    W_dtm = NaN(1,2);
    idc_list = [id1, id2];
    for jj = 1:2
        j = idc_list(jj);
        D_shields_all(j) = D0_S;
        opts.h_bf_all = h_bf_all;
        opts.h_min = (tau_star_critical * delta * d_50) / S0;
        [W_dtm(jj), left_pt_j, right_pt_j, A_bf_all(j), ~, level_code_all(j)] = W_bankfull_from_dtm_modifiche( ...
            j, z_rect_all, D0_S, Z, Xgrid, Ygrid, all_centerlines, opts);
        left_pt_all(j,:) = left_pt_j;
        right_pt_all(j,:) = right_pt_j;
        hydro_data(j).A_bf = A_bf_all(j);
        hydro_data(j).left_pt = left_pt_j;
        hydro_data(j).right_pt = right_pt_j;
    end
    Wm_br = reshape(widths_all([id1 id2]), 1, []);
    if strcmp(node_width_source, 'mask')
        W_dtm = Wm_br;                   % sensitivity test: widths of the flight mask
    end
    W_bf_dtm_all([id1 id2]) = W_dtm;
    if all(isfinite(W_dtm)) && all(W_dtm > 0)
        W1 = W0 * W_dtm(1) / sum(W_dtm);
    elseif all(isfinite(Wm_br)) && all(Wm_br > 0)
        warning('Bif %d: invalid bankfull widths of the branches -> W1/W2 from the mask widths', section_mother);
        W1 = W0 * Wm_br(1) / sum(Wm_br);
    else
        warning('Bif %d: invalid bankfull widths of the branches -> W1 = W2 = W0/2', section_mother);
        W1 = W0 / 2;
    end
    W2 = W0 - W1;
    W_bf_all([id1 id2]) = [W1 W2];
    hydro_data(id1).W_bf = W1;
    hydro_data(id2).W_bf = W2;

    % --- Share of branch 1 in the hydraulic module at flight conditions ---
    % From the depths of the hydraulic module at the first section of each
    % branch, h = z_DTM - z_rect_all, and uniform flow Q ~ W*h^(5/3)*sqrt(S):
    % used only to choose among the asymmetric solutions (see below).
    h_hyd = all_centerlines([id1 id2], 4) - z_rect_all([id1 id2]);
    Qh = reshape(widths_all([id1 id2]), [], 1) .* max(h_hyd, 0).^(5/3) .* sqrt(max(slope_sect([id1 id2]), 0));
    if all(isfinite(Qh)) && sum(Qh) > 0
        psi_hyd = Qh(1) / sum(Qh);
    else
        psi_hyd = NaN;
    end
    psi_hyd_all(i) = psi_hyd;

    % --- Known quantities of the node ---
    q0 = Q_bf / W0;               % bankfull discharge per unit width [m^2/s]
    eta0 = z_rect_all(i);         % bed elevation of the mother (hydraulic module)
    [qs0f_bf, qs0c_bf] = local_qs_fine_coarse(rho, g, S0, D0, ...
        d_f, d_c, d_ref, f_i_f, teta_r, delta);       % bankfull capacity of the mother = sediment supply
    qsf_bf_all(i) = qs0f_bf; qsc_bf_all(i) = qs0c_bf;
    % Gravitational term (§4.4, "grain independence"):
    % theta_x = theta_g0 * l_hr,x with l_hr,x = d_g0/d_x  =>  theta_x = tau_bf0/(rho g delta d_x)
    tau_bf0 = rho * g * D0 * S0;
    d_g0 = local_dg(f_i_f, d_f, d_c, d_ref);
    tetag0 = tau_bf0 / (rho * g * delta * d_g0);
    theta_f = tetag0 * (d_g0 / d_f);
    theta_c = tetag0 * (d_g0 / d_c);

    % --- Node system in the BRT / Ragno form: 6 equations in [psi, D1, D2, S_b, f1F, f2F] ---
    % (see local function morpho_system_brt and the header of the script)
    psi_sol = NaN; Sb = NaN; f1f = NaN; f2f = NaN; D1 = NaN; D2 = NaN;
    inputs_bf = [Q_bf D0 W0 W1 W2 S0 qs0f_bf qs0c_bf f_i_f tetag0];
    if ~all(isfinite(inputs_bf)) || any([Q_bf D0 W0 W1 W2 S0 qs0f_bf qs0c_bf] <= 0)
        warning(['Bif %d: incomplete bankfull data (W0=%.3g, W1=%.3g, W2=%.3g, Q_bf=%.3g) ' ...
                 '-> node system not solved'], section_mother, W0, W1, W2, Q_bf);
    else
        F = @(x) morpho_system_brt(x, Q_bf, D0, W0, W1, W2, S0, ...
                    qs0f_bf, qs0c_bf, f_i_f, rho, g, delta, d_f, d_c, d_ref, teta_r, ...
                    alpha, r, theta_f, theta_c, n_manning);
        % Multi-start: balanced configuration first (psi = W1/W0, D1 = D2 = D0,
        % S_b = S0), then asymmetric ones. The slope is solved as log(S_b).
        X0 = [W1/W0, D0, D0, log(S0), f_i_f, f_i_f];
        for ps = psi_start
            for a1 = D_start
                for a2 = D_start
                    X0(end+1, :) = [ps, a1*D0, a2*D0, log(S0), f_i_f, f_i_f]; %#ok<AGROW>
                end
            end
        end
        sols = zeros(0, 6);                % [psi D1 D2 S_b f1F f2F] of the physical solutions
        all_sols = zeros(0, 9);            % every converged solution, physical or not (diagnostics)
        best_res = Inf; best_exit = NaN;
        for k0 = 1:size(X0, 1)
            [x_k, fval_k, exit_k] = fsolve(F, X0(k0,:), options);
            res_k = norm(fval_k);
            if res_k < best_res
                best_res = res_k; best_exit = exit_k;
            end
            psi_k = x_k(1); D1_k = x_k(2); D2_k = x_k(3); Sb_k = exp(x_k(4));
            physical = exit_k > 0 && res_k < tol_res && ...
                       psi_k > psi_margin && psi_k < 1 - psi_margin && D1_k > 0 && D2_k > 0 && ...
                       all(x_k(5:6) >= 0 & x_k(5:6) <= 1);
            if physical && (isempty(sols) || all(abs(sols(:,1) - psi_k) > 1e-3))
                sols(end+1, :) = [psi_k, D1_k, D2_k, Sb_k, x_k(5), x_k(6)]; %#ok<AGROW>
            end
            % Diagnostics: every converged solution with the reason why it is
            % (not) accepted: 0 = physical, 1 = psi outside (psi_margin, 1 - psi_margin),
            % 2 = D1 or D2 <= 0, 3 = f1F or f2F outside [0, 1].
            if exit_k > 0 && res_k < tol_res && ...
               (isempty(all_sols) || all(abs(all_sols(:,1) - psi_k) > 1e-3))
                if physical
                    why = 0;
                elseif ~(psi_k > psi_margin && psi_k < 1 - psi_margin)
                    why = 1;
                elseif ~(D1_k > 0 && D2_k > 0)
                    why = 2;
                else
                    why = 3;
                end
                all_sols(end+1, :) = [psi_k, D1_k, D2_k, Sb_k / S0, D0 - D1_k, ...
                                      x_k(5), x_k(6), res_k, why]; %#ok<AGROW>
            end
        end
        bif_exitflag_all(i) = best_exit;
        bif_resnorm_all(i) = best_res;
        bif_nsol_all(i) = size(sols, 1);
        node_solutions{i} = sortrows(all_sols, 1);
        if isempty(sols)
            warning('Bif %d: no physical solution of the node system (best residual %.2e, exitflag %d)', ...
                section_mother, best_res, best_exit);
        else
            % Choice of the solution (stability, Bolla Pittaluga et al. 2003):
            % the balanced solution (psi = W1/W0, D1 = D2 = D0, S_b = S0) always
            % exists. Below the critical aspect ratio it is the only one and it
            % is stable. Above it, asymmetric solutions appear and the balanced
            % one becomes unstable: an asymmetric solution is taken, the one
            % closest to the share of branch 1 in the hydraulic module at
            % flight conditions (psi_hyd), i.e. the equilibrium state that is
            % consistent with the branch the data indicate as dominant.
            is_bal = abs(sols(:,1) - W1/W0) < 1e-3 & abs(sols(:,2) - D0) < 1e-3*D0 & ...
                     abs(sols(:,3) - D0) < 1e-3*D0;
            asym = sols(~is_bal, :);
            if isempty(asym)
                if any(is_bal), sel = sols(find(is_bal, 1), :); else, sel = sols(1, :); end
                bif_type_all(i) = 1;
            else
                if isfinite(psi_hyd)
                    [~, k] = min(abs(asym(:,1) - psi_hyd));
                elseif Wm_br(1) >= Wm_br(2)
                    [~, k] = max(asym(:,1));       % no hydraulic share: wider mask branch dominant
                else
                    [~, k] = min(asym(:,1));
                end
                sel = asym(k, :);
                bif_type_all(i) = 2;
            end
            psi_sol = sel(1); D1 = sel(2); D2 = sel(3); Sb = sel(4);
            f1f = sel(5); f2f = sel(6);
            psi_all(i) = psi_sol;
            Sb_all(i) = Sb;
        end
    end

    % Bankfull quantities of the branches (eta_i = eta0 + D0 - D_i: common
    % water surface at the node; stored for diagnostics)
    eta1_sol = eta0 + D0 - D1;
    eta2_sol = eta0 + D0 - D2;
    q1_all(id1) = psi_sol * Q_bf / W1;   q2_all(id2) = (1 - psi_sol) * Q_bf / W2;
    eta1_all(id1) = eta1_sol;  eta2_all(id2) = eta2_sol;
    h1_bf_all(id1) = D1;       h2_bf_all(id2) = D2;
    h_bf_all(id1) = D1;        h_bf_all(id2) = D2;
    Q_bf_all(id1) = psi_sol * Q_bf;
    Q_bf_all(id2) = (1 - psi_sol) * Q_bf;
    hydro_data(id1).h_bf = D1; hydro_data(id2).h_bf = D2;
    hydro_data(id1).Q_bf = Q_bf_all(id1); hydro_data(id2).Q_bf = Q_bf_all(id2);
    f1f_all(id1) = f1f; f1c_all(id1) = 1 - f1f;
    f2f_all(id2) = f2f; f2c_all(id2) = 1 - f2f;

    if isfinite(psi_sol)
        [qs1f, qs1c] = local_qs_fine_coarse(rho, g, Sb, D1, d_f, d_c, d_ref, f1f, teta_r, delta);
        [qs2f, qs2c] = local_qs_fine_coarse(rho, g, Sb, D2, d_f, d_c, d_ref, f2f, teta_r, delta);
        qsf_bf_all([id1 id2]) = [qs1f qs2f];
        qsc_bf_all([id1 id2]) = [qs1c qs2c];
        % bankfull split ratios: water (psi) and sediment per class (psi_x)
        share_Q  = psi_sol;
        share_Qs = [f1f * qs1f * W1 / (f_i_f * qs0f_bf * W0), ...
                    (1 - f1f) * qs1c * W1 / (f_i_c * qs0c_bf * W0)];
        f_seg_f(s1) = f1f;
        f_seg_f(s2) = f2f;
    else
        % Fallback (no physical solution): split proportional to the
        % widths, as in the hydraulic module; branches keep the mother composition.
        Wm = widths_all([id1 id2]);
        if isfinite(W1) && isfinite(W0) && W0 > 0
            share_Q = W1 / W0;
        elseif all(isfinite(Wm)) && all(Wm > 0)
            share_Q = Wm(1) / sum(Wm);
        else
            share_Q = 0.5;
        end
        share_Qs = [share_Q, share_Q];
        f_seg_f([s1 s2]) = f_i_f;
        bif_type_all(i) = 0;
        psi_all(i) = share_Q;
        warning('Bif %d: split proportional to the widths (branch 1 share = %.2f)', ...
            section_mother, share_Q);
    end

    % --- Back to the conditions of the LiDAR flight (§4.3) ---
    % The bankfull split psi is applied to the flight discharge of the mother
    % and psi_x to the sediment fluxes, so water and sediment are conserved:
    %   Q1 = psi*Q0,  Q2 = (1-psi)*Q0,  Q_s1,x = psi_x*Q_s0,x
    % The first section of each branch is then computed in the loop with
    % uniform flow (section B), like every other section.
    Q_map(id1) = Q_map(i) * share_Q;
    Q_map(id2) = Q_map(i) - Q_map(id1);
    q_norm_all(id1) = Q_map(id1);
    q_norm_all(id2) = Q_map(id2);
    Qs_seg(s1,:) = Qs_seg(s_i,:) .* share_Qs;
    Qs_seg(s2,:) = Qs_seg(s_i,:) - Qs_seg(s1,:);
    % Branch widths scaled by W0_mask/W0_bf (diagnostics: the flight
    % hydraulics of the branches uses the mask widths, as downstream).
    ratio_W = W_i / W0;
    for idc = [id1, id2]
        if isfinite(ratio_W)
            W_norm_all(idc) = W_bf_all(idc) * ratio_W;
        else
            W_norm_all(idc) = widths_all(idc);
        end
        hydro_data(idc).W_norm = W_norm_all(idc);
        hydro_data(idc).q_norm = q_norm_all(idc);
    end
end


%% 5. Consistency checks
% 1) Water discharge conservation at the nodes (flight conditions)
for b = 1:size(bif_info,1)
    im = section_to_index(bif_info{b,1});
    ch = section_to_index(bif_info{b,2});
    fprintf('Bifurcation %d: Q mother = %.3f, Q1 = %.3f, Q2 = %.3f m^3/s (error = %.1e)\n', ...
        bif_info{b,1}, Q_map(im), Q_map(ch(1)), Q_map(ch(2)), Q_map(im) - Q_map(ch(1)) - Q_map(ch(2)));
end
for c = 1:size(conf_info,1)
    ic = section_to_index(conf_info{c,1});
    up = section_to_index(conf_info{c,2});
    fprintf('Confluence %d: Q downstream = %.3f, Q1 = %.3f, Q2 = %.3f m^3/s\n', ...
        conf_info{c,1}, Q_map(ic), Q_map(up(1)), Q_map(up(2)));
end

% 2) Sediment flux conservation per class at the bifurcations, and node solution
for b = 1:size(bif_info,1)
    im = section_to_index(bif_info{b,1});
    ch = section_to_index(bif_info{b,2});
    fprintf('Bifurcation %d: Qs,F mother = %.3e, branches = %.3e + %.3e | Qs,C mother = %.3e, branches = %.3e + %.3e m^3/s\n', ...
        bif_info{b,1}, qs_tot_f_all(im), qs_tot_f_all(ch(1)), qs_tot_f_all(ch(2)), ...
        qs_tot_c_all(im), qs_tot_c_all(ch(1)), qs_tot_c_all(ch(2)));
    fprintf(['   bankfull: W0 = %.2f, W1 = %.2f, W2 = %.2f m | q1 = %.3f, q2 = %.3f m^2/s | D1 = %.2f, D2 = %.2f m | ' ...
             'S_b/S0 = %.2f | eta1-eta0 = %.3f, eta2-eta0 = %.3f m | f1F = %.3f, f2F = %.3f\n'], ...
        W_bf_all(im), W_bf_all(ch(1)), W_bf_all(ch(2)), q1_all(ch(1)), q2_all(ch(2)), ...
        h1_bf_all(ch(1)), h2_bf_all(ch(2)), Sb_all(im) / S_node_all(im), ...
        eta1_all(ch(1)) - z_rect_all(im), eta2_all(ch(2)) - z_rect_all(im), f1f_all(ch(1)), f2f_all(ch(2)));
end

% 3) Diagnostics of the bifurcation nodes, in three tables.
%  a) Slopes and depths. S0: slope of the mother used in the node (window
%     upstream of the node); S1m, S2m: slopes measured downstream of the node
%     (comparison only, they are not used); Sb: equilibrium slope of the
%     branches from the node system; src: origin of S0, S1m, S2m (1 = window
%     L, 2 = enlarged window, 3 = fallback, 4 = raised to the floor);
%     S0seg: segment slope slope_sect of the mother; D0S: Shields depth of
%     the mother; D0: depth used; lvl: 0 = Shields, 1 = limited by the banks,
%     2 = spill level below the critical depth -> D0 = critical depth (low
%     banks); th/thc: bankfull Shields number of the mother over the critical one.
fprintf('\nNode slopes and bankfull depths (node_width_source = %s)\n', node_width_source);
fprintf('%-6s %-9s %-9s %-9s %-9s %-5s %-9s %-7s %-7s %-5s %-6s\n', ...
    'Bif', 'S0', 'S1m', 'S2m', 'Sb', 'src', 'S0seg', 'D0S', 'D0', 'lvl', 'th/thc');
for b = 1:size(bif_info,1)
    im = section_to_index(bif_info{b,1});
    ch = section_to_index(bif_info{b,2});
    S_b = S_node_all([im ch(1) ch(2)]);
    src_b = S_node_src_all([im ch(1) ch(2)]);
    fprintf('%-6d %-9.2e %-9.2e %-9.2e %-9.2e %-5s %-9.2e %-7.2f %-7.2f %-5d %-6.2f\n', ...
        bif_info{b,1}, S_b(1), S_b(2), S_b(3), Sb_all(im), sprintf('%d%d%d', src_b), ...
        slope_sect(im), D_shields_all(im), h_bf_all(im), level_code_all(im), ...
        theta_bf_all(im) / tau_star_critical);
end
%  b) Widths and node solution. W0m, W1m, W2m: mask widths; W0bf: width of
%     the mother used in the node; W1dtm, W2dtm: branch widths before the
%     rescaling on W0; beta = W0/(2*D0); Qbf: bankfull discharge of the
%     mother; psi: selected share of the discharge to branch 1; W1/W0: share
%     of the balanced solution; psi_hyd: share of branch 1 in the hydraulic
%     module at flight conditions (used to choose among asymmetric
%     solutions); nsol: number of distinct physical solutions; type:
%     fallback / balanced / asymmetric.
type_names = {'fallback', 'balanced', 'asymmetric'};
fprintf('\nNode widths and solution\n');
fprintf('%-6s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-8s %-6s %-6s %-7s %-5s %-9s %s\n', ...
    'Bif', 'W0m', 'W0bf', 'W1m', 'W2m', 'W1dtm', 'W2dtm', 'beta', 'Qbf', 'psi', 'W1/W0', ...
    'psi_hyd', 'nsol', 'resnorm', 'type');
for b = 1:size(bif_info,1)
    im = section_to_index(bif_info{b,1});
    ch = section_to_index(bif_info{b,2});
    if isfinite(bif_type_all(im)), tname = type_names{bif_type_all(im) + 1}; else, tname = 'not computed'; end
    fprintf('%-6d %-7.1f %-7.1f %-7.1f %-7.1f %-7.1f %-7.1f %-7.1f %-8.1f %-6.2f %-6.2f %-7.2f %-5g %-9.2e %s\n', ...
        bif_info{b,1}, widths_all(im), W_bf_all(im), widths_all(ch(1)), widths_all(ch(2)), ...
        W_bf_dtm_all(ch(1)), W_bf_dtm_all(ch(2)), beta_all(im), Q_bf_all(im), psi_all(im), ...
        W_bf_all(ch(1)) / W_bf_all(im), psi_hyd_all(im), bif_nsol_all(im), bif_resnorm_all(im), tname);
end
%  c) All converged solutions of each node system (physical or not).
%     psi = share of the discharge to branch 1; D1, D2 = bankfull depths of
%     the branches; Sb/S0 = equilibrium slope of the branches over the slope
%     of the mother; deta = eta1 - eta0 = D0 - D1; status: ok, or the reason
%     why the solution is discarded.
status_names = {'ok', 'psi out of range', 'depth <= 0', 'f out of [0,1]'};
fprintf('\nAll converged solutions of the node systems\n');
for b = 1:size(bif_info,1)
    im = section_to_index(bif_info{b,1});
    ch = section_to_index(bif_info{b,2});
    sl = node_solutions{im};
    fprintf('Bif %d: W1/W0 = %.2f, psi_hyd = %.2f, beta = %.1f, selected psi = %.2f (%d converged solutions)\n', ...
        bif_info{b,1}, W_bf_all(ch(1)) / W_bf_all(im), psi_hyd_all(im), beta_all(im), psi_all(im), size(sl, 1));
    if isempty(sl)
        continue
    end
    fprintf('   %-7s %-7s %-7s %-7s %-8s %-6s %-6s %s\n', ...
        'psi', 'D1', 'D2', 'Sb/S0', 'deta', 'f1F', 'f2F', 'status');
    for k = 1:size(sl, 1)
        fprintf('   %-7.3f %-7.2f %-7.2f %-7.3f %-8.3f %-6.3f %-6.3f %s\n', ...
            sl(k,1), sl(k,2), sl(k,3), sl(k,4), sl(k,5), sl(k,6), sl(k,7), status_names{sl(k,9) + 1});
    end
end

% 4) Sections without a morphodynamic bed elevation
no_bed = find(isnan(z_rect_morpho));
fprintf('%d of %d sections without a morphodynamic bed elevation\n', numel(no_bed), N);
if ~isempty(no_bed)
    fprintf('  segment codes involved: %s\n', mat2str(unique(all_centerlines(no_bed,6))'));
end

% 5) Comparison with the hydraulic module
dz_morpho = z_rect_morpho - z_rect_all;
fprintf('z_rect_morpho - z_rect_all: mean = %.3f m, max |dz| = %.3f m\n', ...
    mean(dz_morpho, 'omitnan'), max(abs(dz_morpho), [], 'omitnan'));


%% 6. Save the results of the morphodynamic module
out_dir = '\\nas-camporeale\home\ELISA\BATHIMETRY\0_MODIFICHE';

% Centerline with the morphodynamic bed elevation (same CSV format as the hydraulic module)
fileID = fopen(fullfile(out_dir, '0_centerline_bat_rect_2025_morpho_BRT.csv'), 'w');
fprintf(fileID, 'X,Y,Z\n');
for i = 1:N
    fprintf(fileID, '%.3f,%.3f,%.3f\n', all_centerlines(i,2), all_centerlines(i,3), z_rect_morpho(i));
end
fclose(fileID);

% Depth and bed elevation of each cross-section
save(fullfile(out_dir, 'h_rect_all_morpho_BRT.mat'), 'h_rect_morpho');
save(fullfile(out_dir, 'z_rect_all_morpho_BRT.mat'), 'z_rect_morpho');

% Bed elevation assigned to every point of each cross-section (same layout as z_matrix)
z_rect_morpho_matrix = NaN(size(z_matrix));
for i = 1:size(z_matrix, 1)
    pts = ~isnan(z_matrix(i, :));
    z_rect_morpho_matrix(i, pts) = z_rect_morpho(i);
end

% All results of the morphodynamic module
save(fullfile(out_dir, 'results_2025_morpho_BRT.mat'), ...
    'z_rect_morpho', 'h_rect_morpho', 'z_rect_morpho_matrix', 'Q_map', ...
    'W_bf_all', 'W_bf_dtm_all', 'h_bf_all', 'Q_bf_all', 'A_bf_all', 'left_pt_all', 'right_pt_all', ...
    'W_norm_all', 'q_norm_all', 'q1_all', 'q2_all', 'h1_bf_all', 'h2_bf_all', 'eta1_all', 'eta2_all', ...
    'f1f_all', 'f1c_all', 'f2f_all', 'f2c_all', 'bif_exitflag_all', 'bif_resnorm_all', 'S_node_all', 'S_node_src_all', 'Sb_all', 'psi_hyd_all', 'node_width_source', 'D_shields_all', 'level_code_all', 'node_solutions', ...
    'beta_all', 'theta_bf_all', 'psi_all', 'bif_type_all', 'bif_nsol_all', ...
    'qs0f_all', 'qs0c_all', 'qsf_cap_all', 'qsc_cap_all', 'qsf_bf_all', 'qsc_bf_all', ...
    'qs_tot_f_all', 'qs_tot_c_all', 'total_QS_all_f', 'total_QS_all_c', 'f0f_all', 'f0c_all', ...
    'hydro_data', 'sediment_data', 'sectionWet');
disp('Morphodynamic results saved.');


%% 7. Bathymetric correction of the DTM with the morphodynamic bed (same procedure as the hydraulic module)
% The bed elevation of each section is assigned to all points of its
% cross-section, interpolated on the DTM grid (scattered interpolant) and
% written inside the wet-area mask of the raw LiDAR DTM, never above the
% original LiDAR elevation. Outputs: riverbed elevation along the
% cross-sections, corrected DTM and correction map (DoD), all with the _morpho
% suffix.
%% ======================= DTM CORRECTION: INPUT RASTERS =======================
% Wet-area mask of the LiDAR flight and RAW LiDAR DTM. The morphodynamic bed
% replaces the raw DTM inside the wet area (same procedure as the hydraulic
% module), so the correction starts from the raw DTM, not from the
% hydraulically corrected one loaded above.
[mask_wetarea, Rmask_wetarea] = readgeoraster('I_wetarea_2025.tif');
mask_wetarea = logical(mask_wetarea);
[raster_dtm, R_dtm] = readgeoraster('cutDTM_50cm_ORCO_UTM32N_ETRF2000.tif');
raster_dtm = double(raster_dtm);
raster_dtm(raster_dtm == -32767) = NaN;   % NoData of the raw DTM

% Are the DTM and the mask aligned?
isequal(R_dtm, Rmask_wetarea)
fprintf('--- DTM ---\n');
disp(R_dtm)
fprintf('--- Wet Area Mask ---\n');
disp(Rmask_wetarea)

% DTM bounding box
xmin_dtm = R_dtm.XWorldLimits(1);
xmax_dtm = R_dtm.XWorldLimits(2);
ymin_dtm = R_dtm.YWorldLimits(1);
ymax_dtm = R_dtm.YWorldLimits(2);
% Mask bounding box
xmin_mask = Rmask_wetarea.XWorldLimits(1);
xmax_mask = Rmask_wetarea.XWorldLimits(2);
ymin_mask = Rmask_wetarea.YWorldLimits(1);
ymax_mask = Rmask_wetarea.YWorldLimits(2);

fprintf('DTM extent: xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f\n', xmin_dtm, xmax_dtm, ymin_dtm, ymax_dtm);
fprintf('Mask extent: xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f\n', xmin_mask, xmax_mask, ymin_mask, ymax_mask);

%% Align the wet-area mask with the DTM grid
% World coordinates of every mask pixel
[cols_mask, rows_mask] = meshgrid(1:size(mask_wetarea,2), 1:size(mask_wetarea,1));
[xw_mask, yw_mask] = Rmask_wetarea.intrinsicToWorld(cols_mask, rows_mask);

% Keep only the TRUE pixels of the mask
mask_true_idx = find(mask_wetarea);
x_mask_true = xw_mask(mask_true_idx);
y_mask_true = yw_mask(mask_true_idx);

% Convert their world coordinates into intrinsic (column, row) indices of the DTM
[col_in_dtm, row_in_dtm] = R_dtm.worldToIntrinsic(x_mask_true, y_mask_true);
col_in_dtm = round(col_in_dtm);
row_in_dtm = round(row_in_dtm);

% Discard indices outside the DTM
valid = col_in_dtm >= 1 & col_in_dtm <= size(raster_dtm,2) & ...
        row_in_dtm >= 1 & row_in_dtm <= size(raster_dtm,1);
col_in_dtm = col_in_dtm(valid);
row_in_dtm = row_in_dtm(valid);

% Wet-area mask on the DTM grid
mask_extended = false(size(raster_dtm));
lin_idx = sub2ind(size(raster_dtm), row_in_dtm, col_in_dtm);
mask_extended(lin_idx) = true;

% Check: sample some TRUE pixels and verify that the coordinates match after
% the round trip mask -> DTM indices -> world coordinates
sample_idx = mask_true_idx(1:100:end);  % one pixel every 100
x_check = xw_mask(sample_idx);
y_check = yw_mask(sample_idx);

[col_check, row_check] = R_dtm.worldToIntrinsic(x_check, y_check);
[x_back, y_back] = R_dtm.intrinsicToWorld(round(col_check), round(row_check));

max_err_x = max(abs(x_check - x_back));
max_err_y = max(abs(y_check - y_back));
fprintf('Max roundtrip error: X=%.4f m, Y=%.4f m\n', max_err_x, max_err_y);

%% Riverbed elevation along the cross-sections (morphodynamic module)
mask_valid = ~isnan(z_rect_morpho_matrix);
x = x_matrix(mask_valid);
y = y_matrix(mask_valid);
z = z_rect_morpho_matrix(mask_valid);
figure('Name', 'Riverbed_elevation_morpho_BRT');
scatter(x, y, 10, z, 'filled');
clim_riverbed = prctile(z, [1 99]);   % percentiles instead of min/max: robust to isolated outliers
caxis(clim_riverbed);
ylabel(colorbar, 'Z (m)');
xlabel('X'); ylabel('Y');
title('Riverbed elevation along cross-sections_morpho_BRT', 'Interpreter', 'none');
axis equal;

%% Interpolate the bed elevation on the DTM grid and correct the DTM
F_bed = scatteredInterpolant(x, y, z, 'natural', 'none');

% Interpolation grid = DTM cell centres
[cols, rows] = meshgrid(1:size(raster_dtm,2), 1:size(raster_dtm,1));
[xq, yq] = R_dtm.intrinsicToWorld(cols, rows);
zq = F_bed(xq, yq);

% Correction inside the wet area only
raster_dtm_corrected_morpho = raster_dtm;

% Never raise the bed above the original LiDAR elevation: the model only adds
% bathymetric information (deepening), it never overwrites upwards what the
% LiDAR has already measured.
n_would_raise = nnz(mask_extended & (zq > raster_dtm));
fprintf('Pixels where the correction would raise the bed above the raw DTM: %d (%.2f%% of the corrected area)\n', ...
    n_would_raise, 100*n_would_raise/nnz(mask_extended));

zq_safe = min(zq, raster_dtm);
raster_dtm_corrected_morpho(mask_extended) = zq_safe(mask_extended);

%% DTM after the morphodynamic bathymetric correction
% Each row and column is placed at its real world coordinate (first and
% last row/column from intrinsicToWorld) and the axis is set to 'xy': north
% is up with correct tick labels, whatever the row order of the raster
% (north-first or south-first).
[x_img, ~] = R_dtm.intrinsicToWorld([1 size(raster_dtm,2)], [1 1]);
[~, y_img] = R_dtm.intrinsicToWorld([1 1], [1 size(raster_dtm,1)]);
figure('Name', 'DTM_corrected_morpho_BRT');
imAlpha = true(size(raster_dtm_corrected_morpho));
imAlpha(raster_dtm_corrected_morpho == 0 | isnan(raster_dtm_corrected_morpho)) = false;
imagesc(x_img, y_img, raster_dtm_corrected_morpho, 'AlphaData', double(imAlpha));
axis xy; axis image;
colormap(parula);
valid_vals = raster_dtm_corrected_morpho(~isnan(raster_dtm_corrected_morpho));
caxis(prctile(valid_vals, [1 99]));
cb = colorbar;
ylabel(cb, 'Z (m)');
xlabel('X'); ylabel('Y');
title('DTM after bathymetric correction_morpho_BRT', 'Interpreter', 'none');

%% Bathymetric correction map (DoD = raw DTM - corrected DTM)
DoD_morpho = raster_dtm - raster_dtm_corrected_morpho;

figure('Name', 'Bathymetric_correction_map_morpho_BRT');
imAlpha = true(size(DoD_morpho));
imAlpha(DoD_morpho == 0 | isnan(DoD_morpho)) = false;   % transparent where no correction was applied
% x_img, y_img (world coordinates of the first/last column and row) come from the previous section
imagesc(x_img, y_img, DoD_morpho, 'AlphaData', double(imAlpha));
axis xy; axis image;
set(gca, 'Color', 'white');           % white background under the transparent pixels
colormap(parula);
nz = DoD_morpho(DoD_morpho > 0 & ~isnan(DoD_morpho));
L = prctile(nz, 98);                  % dynamic limit, robust to isolated extreme pixels
caxis([0 L]);
cb = colorbar;
ylabel(cb, 'ΔZ (m)');
title('Bathymetric correction map_morpho_BRT', 'Interpreter', 'none');
xlabel('X'); ylabel('Y');

% Save the corrected DTM and the DoD as GeoTIFF (EPSG:32632)
geotiffwrite(fullfile(out_dir, '0_DTM_bat_corrected_2025_morpho_BRT.tif'), raster_dtm_corrected_morpho, R_dtm, 'CoordRefSysCode', 32632);
geotiffwrite(fullfile(out_dir, '0_DoD_bat_corrected_2025_morpho_BRT.tif'), DoD_morpho, R_dtm, 'CoordRefSysCode', 32632);
disp('Done!');


%% 8. Local functions
% =====================================================================
%  Local functions (local functions in a script must be at the end of the file)
%  External functions on the path: W_bankfull_from_dtm_modifiche.m, G_zi.m
% =====================================================================

% ==============  Bankfull node system, BRT / Ragno form  ==============
function R = morpho_system_brt(x, Q0, D0, W0, W1, W2, S0, ...
                qs0f, qs0c, f0f, rho, g, delta, d_f, d_c, d_ref, teta_r, ...
                alpha, r, theta_f, theta_c, n_manning)
% Unknowns x = [psi, D1, D2, log(S_b), f1F, f2F]:
%   psi  = Q1/Q0, share of the bankfull discharge to branch 1 (Q2 = Q0 - Q1)
%   D1, D2 = bankfull depths of the branches
%   S_b  = equilibrium slope of the branches (common to both branches: equal
%          water level at the node and at the downstream end of branches of
%          equal length, Bolla Pittaluga et al. 2003). Solved as log(S_b) > 0.
%   f1F, f2F = fine fraction of the bed surface of the branches
% Equations (6):
%   (1)-(2) uniform flow in each branch:  Q_i = W_i * D_i^(5/3) * sqrt(S_b) / n
%   (3)-(4) cell balance for F and C (two cells of length alpha*W0 upstream
%           of the node), with the transverse bed slope (eta1 - eta2)/(W0/2)
%           and eta1 - eta2 = D2 - D1 (common water surface at the node):
%             q_sy,x = q_s0,x * ( q_y/q0 - r/sqrt(theta_x) * (eta1 - eta2)/(W0/2) )
%             q_y    = Q_y/(alpha*W0),  Q_y = (Q1 - Q2 - Q0*(W1 - W2)/W0)/2
%             f0x*q_sy,x * alpha*W0 = ( Qs1x - Qs2x - Qs0x*(W1-W2)/(W1+W2) )/2
%   (5)-(6) sediment continuity for F and C: Qs1x + Qs2x = Qs0x
% Water continuity Q1 + Q2 = Q0 is built in. The bed elevation of the
% node, eta0 = (eta1 + eta2)/2 of the previous formulation, is NOT used:
% the equilibrium slope of the branches is the unknown that closes the system.
% qs0f, qs0c = bankfull transport capacity per unit width of the mother
% (sediment supply); theta_x = theta_g0 * d_g0/d_x (bankfull Shields number
% of class x in the mother).

    psi = x(1);
    D1 = max(x(2), 1e-6);
    D2 = max(x(3), 1e-6);
    Sb = exp(x(4));
    f1F = x(5);
    f2F = x(6);
    Q1 = psi * Q0;
    Q2 = Q0 - Q1;

    % (1)-(2) uniform flow, non-dimensional
    uf1 = Q1 / (W1 * D1^(5/3) * sqrt(Sb) / n_manning) - 1;
    uf2 = Q2 / (W2 * D2^(5/3) * sqrt(Sb) / n_manning) - 1;

    % transport capacity per unit width (Wilcock-Crowe, d_g of each branch)
    [qs1f, qs1c] = local_qs_fine_coarse(rho, g, Sb, D1, d_f, d_c, d_ref, f1F, teta_r, delta);
    [qs2f, qs2c] = local_qs_fine_coarse(rho, g, Sb, D2, d_f, d_c, d_ref, f2F, teta_r, delta);

    f0c = 1 - f0f;
    QS0f = f0f * qs0f * W0;          QS0c = f0c * qs0c * W0;
    QS1f = f1F * qs1f * W1;          QS1c = (1 - f1F) * qs1c * W1;
    QS2f = f2F * qs2f * W2;          QS2c = (1 - f2F) * qs2c * W2;

    % (5)-(6) sediment continuity, non-dimensional
    cons_f = (QS1f + QS2f) / QS0f - 1;
    cons_c = (QS1c + QS2c) / QS0c - 1;

    % (3)-(4) cell balance, divided by f0x*q_s0x*alpha (units of m)
    q0 = Q0 / W0;
    Qy_over_q0 = (Q1 - Q2 - Q0*(W1 - W2)/W0) / (2*q0);       % Q_y / q0  [m]
    slope_y = (D2 - D1) / (W0/2);                              % (eta1 - eta2)/(W0/2)
    eq_f = Qy_over_q0/alpha - W0*r/sqrt(theta_f)*slope_y ...
           - (QS1f - QS2f - QS0f*(W1 - W2)/(W1 + W2)) / (2*alpha*f0f*qs0f);
    eq_c = Qy_over_q0/alpha - W0*r/sqrt(theta_c)*slope_y ...
           - (QS1c - QS2c - QS0c*(W1 - W2)/(W1 + W2)) / (2*alpha*f0c*qs0c);
    % scaled by W0 so that all residuals are O(1)
    R = [uf1; uf2; eq_f / W0; eq_c / W0; cons_f; cons_c];
end

function m = cell2mat_safe(c)
% First column of a cell array as a numeric vector ([] if the cell is empty)
    if isempty(c)
        m = [];
    else
        m = cell2mat(c(:,1));
    end
end

function d_g = local_dg(f_F, d_f, d_c, d_ref)
% Geometric mean diameter of the bed surface (methodology §4.1):
%   phi_m = f_C*phi_C + f_F*phi_F,  d_g = d_ref * 2^(-phi_m)
    f_F = min(max(f_F, 0), 1);
    phi_f = -log2(d_f / d_ref);
    phi_c = -log2(d_c / d_ref);
    phi_m = (1 - f_F) * phi_c + f_F * phi_f;
    d_g = d_ref * 2^(-phi_m);
end

function [qs_f, qs_c] = local_qs_fine_coarse(rho, g, S, h, d_f, d_c, d_ref, f_F, teta_r, delta)
% Transport capacity per unit width of the fine and coarse classes
% (Wilcock & Crowe 2003, methodology §4.2), to be multiplied by f_x:
%   theta_g = tau/(rho*g*delta*d_g),  zeta_x = g_hr,x * theta_g/theta_r
%   g_hr,x  = (d_x/d_g)^(-b_x),        b_x = 0.67/(1 + exp(1.5 - d_x/d_g))
%   qs_x    = (tau/rho)^(3/2)/(g*delta) * G(zeta_x)
% d_g depends on the surface composition f_F.
    tau_loc = rho * g * S * h;
    d_g = local_dg(f_F, d_f, d_c, d_ref);
    teta_g = tau_loc / (rho * g * delta * d_g);

    b_f = 0.67 / (1 + exp(1.5 - d_f / d_g));
    zi_f = (d_f / d_g)^(-b_f) * teta_g / teta_r;
    qs_f = (tau_loc / rho)^(3/2) / (g * delta) * G_zi(zi_f);

    b_c = 0.67 / (1 + exp(1.5 - d_c / d_g));
    zi_c = (d_c / d_g)^(-b_c) * teta_g / teta_r;
    qs_c = (tau_loc / rho)^(3/2) / (g * delta) * G_zi(zi_c);
end

function f = local_equilibrium_f(QsF, QsC, rho, g, S, h, d_f, d_c, d_ref, teta_r, delta, f_default)
% Surface composition that carries fine and coarse sediment in the ratio of
% the supply (methodology §5.4): f*qsF(f) = R*(1-f)*qsC(f), R = QsF/QsC.
% (f = 0 -> residual < 0, f = 1 -> residual > 0: a root in [0,1] exists)
    f = f_default;
    if ~(isfinite(QsF) && isfinite(QsC) && QsF >= 0 && QsC > 0 && h > 0 && S > 0)
        return
    end
    R = QsF / QsC;
    fun = @(ff) residual_eq_f(ff, R, rho, g, S, h, d_f, d_c, d_ref, teta_r, delta);
    if fun(0) * fun(1) < 0
        f = fzero(fun, [0 1]);
    elseif R == 0
        f = 0;
    end
end

function res = residual_eq_f(ff, R, rho, g, S, h, d_f, d_c, d_ref, teta_r, delta)
    [qf, qc] = local_qs_fine_coarse(rho, g, S, h, d_f, d_c, d_ref, ff, teta_r, delta);
    res = ff * qf - R * (1 - ff) * qc;
end

function ks = local_ks_i(i, ks_cd50, ks_default)
% Strickler coefficient [m^(1/3)/s] of section i:
% empty ks_cd50 -> ks_default (= 1/n_manning); scalar -> that value;
% vector -> ks_cd50(i) if finite, otherwise ks_default.
    if isempty(ks_cd50)
        ks = ks_default;
    elseif isscalar(ks_cd50)
        ks = ks_cd50;
    elseif i <= numel(ks_cd50) && isfinite(ks_cd50(i)) && ks_cd50(i) > 0
        ks = ks_cd50(i);
    else
        ks = ks_default;
    end
end

function h = local_normal_h(Q, i, widths_all, slope_sect, all_centerlines, h_rect_morpho, ks_cd50, ks_default)
% Uniform-flow depth in a wide rectangular section (Gauckler-Strickler):
%   Q = ks * W * h^(5/3) * sqrt(S)   ->   h = (Q / (ks * W * sqrt(S)))^(3/5)
% With ks = 1/n_manning it is the same formula as the hydraulic module.
% If W or S of the section are not valid, the previous section of the same
% segment is used (its geometry, or as a last resort its depth).
    h = NaN;
    if ~isfinite(Q) || Q <= 0
        return
    end
    W = widths_all(i);
    S = slope_sect(i);
    ks = local_ks_i(i, ks_cd50, ks_default);
    if ~(W > 0 && S > 0)
        same_seg_prev = i > 1 && all_centerlines(i-1,6) == all_centerlines(i,6);
        if ~same_seg_prev
            return
        end
        W = widths_all(i-1);
        S = slope_sect(i-1);
        ks = local_ks_i(i-1, ks_cd50, ks_default);
        if ~(W > 0 && S > 0)
            h = h_rect_morpho(i-1);
            return
        end
    end
    h = (Q / (ks * W * sqrt(S)))^(3/5);
end
% ==============  Slopes at the bifurcation nodes (methodology §7.2)  ==============
function [S, src] = local_node_slope(i_start, direction, L, L_min, L_max, n_min, S_floor, S_fallback, net)
% Water-surface slope over a window of length L starting at the node.
%   direction = -1 : window UPSTREAM of i_start (mother, i_start = last
%                    point of the mother segment)
%   direction = +1 : window DOWNSTREAM of i_start (branch, i_start = first
%                    point of the branch segment)
% The window follows the centerline point by point and, when it reaches
% the end of the current segment, continues into the adjacent segment:
%   - upstream: the predecessor; at a confluence (two predecessors) the one
%     with the larger discharge Q_map at its last point (main channel);
%   - downstream: the successor; at a further bifurcation (two successors)
%     the one with the larger mask width at its first point.
% The window stops at the network head / outlet or after length L.
% The slope is S = -dz/ds of the least-squares line z(s) through the
% water-surface elevations all_centerlines(:,4) of the window points
% (s = along-channel distance, positive downstream).
% If the fit is not usable (fewer than n_min points, window shorter than
% L_min because the network ends, or S <= 0) the window is doubled up to
% L_max. If it is still not usable, S_fallback is used (a typical slope of
% the network, 25th percentile of the segment slopes): the segment slope
% slope_sect is not used because at these nodes it is often the 1.6e-4
% floor of the pre-processing.
% A valid window slope is kept as measured; it is only raised to S_floor
% (5th percentile of the segment slopes) to discard values close to zero.
% src: 1 = window L, 2 = enlarged window, 3 = S_fallback (no valid window),
%      4 = window slope raised to S_floor.
    S = NaN;
    src = NaN;
    L_try = L;
    while L_try <= L_max
        [rows, s_dist] = local_node_window(i_start, direction, L_try, net);
        z = net.all_centerlines(rows, 4);
        ok = isfinite(z);
        if nnz(ok) >= n_min && (max(s_dist(ok)) - min(s_dist(ok))) >= L_min
            p = polyfit(s_dist(ok), z(ok), 1);
            if -p(1) > 0
                S = -p(1);
                src = 1 + (L_try > L);
                break
            end
        end
        L_try = 2 * L_try;
    end
    if ~isfinite(S)
        S = S_fallback;
        src = 3;
    elseif S < S_floor
        S = S_floor;
        src = 4;
    end
end

function [rows, s_dist] = local_node_window(i_start, direction, L, net)
% Centerline rows of the window of length L from i_start (see
% local_node_slope) and their along-channel distance s, positive downstream
% and measured from i_start.
    xy = net.all_centerlines(:, 2:3);
    rows = i_start;
    s_dist = 0;
    r = i_start;
    dist = 0;
    for it = 1:size(xy, 1)                       % safety bound on the steps
        seg = net.seg_of_row(r);
        if direction < 0
            if r > net.sections_ids(seg, 1)
                r_next = r - 1;                  % previous point, same segment
            else
                up = net.predecessor{seg};
                if isempty(up)
                    break                        % network head
                elseif numel(up) > 1             % confluence: main (larger discharge) upstream branch
                    [~, k] = max(net.Q_map(net.sections_ids(up, 2)));
                    up = up(k);
                end
                r_next = net.sections_ids(up, 2);
            end
        else
            if r < net.sections_ids(seg, 2)
                r_next = r + 1;                  % next point, same segment
            else
                dn = net.successor{seg};
                if isempty(dn)
                    break                        % network outlet
                elseif numel(dn) > 1             % further bifurcation: wider downstream branch
                    [~, k] = max(net.widths_all(net.sections_ids(dn, 1)));
                    dn = dn(k);
                end
                r_next = net.sections_ids(dn, 1);
            end
        end
        dist = dist + norm(xy(r_next, :) - xy(r, :));
        if dist > L
            break
        end
        r = r_next;
        rows(end+1, 1) = r;                      %#ok<AGROW>
        s_dist(end+1, 1) = direction * dist;     %#ok<AGROW>
    end
end
