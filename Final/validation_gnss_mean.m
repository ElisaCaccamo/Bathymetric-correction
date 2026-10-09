function [Tstats, Tpts] = validation_gnss_mean(gnss_path, raster_dtm, raster_dtm_corrected, R_dtm, mask_extended, out_dir)
% VALIDATION_GNSS_MEAN  Validate the original DTM and the bathymetry-corrected
% DTM against GNSS field points, using the section-mean bed elevation as
% the reference "truth" value for every individual point.
%
% WHAT THIS COMPARES (read this before interpreting the output)
%   1) GNSS points are grouped by measured cross-section (parsed from point names
%      following the pattern "sNNN_pMMM", e.g. "s006_p015" -> station
%      "s006"). For each station, the MEAN bed elevation across its points
%      is computed (X,Y are NOT averaged).
%   2) That single mean elevation is assigned back to EVERY individual
%      point of the station, at its own, unchanged (X,Y) location -
%      stations are not collapsed to a centroid, full point density and
%      spatial resolution are kept.
%   3) The DTM (original and corrected) is sampled, bilinearly, at each
%      point's own real (X,Y) - one sample per point, as in the point-wise
%      validation.
%   4) The residual at each point is:
%         e = DTM_sampled_at_that_point - cross_section_mean_elevation
%      i.e. every point is compared against the mean of its own section,
%      not against its own individual GNSS reading.
%   Rationale: the hydraulic model assigns a single flat elevation to the
%   whole cross-section (z_rect). Comparing that flat value against single,
%   laterally-scattered GNSS readings penalises the model for transverse
%   variability it was never meant to reproduce. Comparing it against the
%   cross-section MEAN removes that mismatch while still keeping one residual per
%   measured point (not per station), preserving statistical power and the
%   full spatial resolution of the residual map.
%
% Residual convention:  e = DTM_sampled_at_that_point - cross_section_mean_elevation
%   e > 0  ->  model ABOVE the mean GNSS elevation
%   e < 0  ->  model BELOW the mean GNSS elevation
%
% INPUTS
%   gnss_path              csv with columns: "NomePunto, Y, X, Quota ellissoidica, Quota ortometrica"
%   raster_dtm             original DTM (double, NaN at nodata)
%   raster_dtm_corrected   bathymetry-corrected DTM (same grid)
%   R_dtm                  DTM spatial reference (MapCellsReference)
%   mask_extended          logical mask (DTM grid) where the correction was applied
%   out_dir                (optional) output folder for csv/figure
%
% OUTPUTS
%   Tstats  summary metrics table (N, ME, Median, SD, RMSE, MAE, NMAD, P95_abs)
%           N = number of POINTS (not cross-sections)
%   Tpts    per-point table: each point keeps its own (X,Y) and DTM sample,
%           but is compared against its cross-section's mean GNSS elevation
%
% Requires MATLAB R2020a+ (VariableNamingRule, tiledlayout, exportgraphics).
% No toolboxes beyond base MATLAB.

if nargin < 6 || isempty(out_dir), out_dir = fileparts(gnss_path); end
rng(1);   % reproducible bootstrap

%% 1) Robust csv read (separator/decimal mark not guaranteed)
opts = detectImportOptions(gnss_path, 'VariableNamingRule', 'preserve');
T = readtable(gnss_path, opts);
if width(T) < 5
    error('GNSS file: expected >= 5 columns, found %d. Check the delimiter.', width(T));
end

Vraw = nan(height(T), 4);
for k = 1:4
    c = T{:, k+1};
    if isnumeric(c)
        Vraw(:, k) = double(c);
    else
        Vraw(:, k) = str2double(strrep(string(c), ',', '.'));   % handle comma decimal mark
    end
end
name  = string(T{:, 1});
Yg    = Vraw(:, 1);   % NB: Y precedes X in the source file
Xg    = Vraw(:, 2);
z_ell = Vraw(:, 3);   % each point's OWN elevation (kept for diagnostics, not used as the comparison truth)
z_ort = Vraw(:, 4);

% Safety check: in UTM 32N (Piedmont) Easting ~4e5 < Northing ~5e6
if median(Xg, 'omitnan') > median(Yg, 'omitnan')
    warning('X and Y look swapped relative to UTM32N expectations: swapping them back.');
    [Xg, Yg] = deal(Yg, Xg);
end

%% 2) Station grouping and per-station MEAN ELEVATION (X,Y are not touched)
station = regexp(name, '^(s\d+)_p\d+$', 'tokens', 'once', 'ignorecase');
has_station = ~cellfun(@isempty, station);
n_missing = nnz(~has_station);
if n_missing > 0
    fprintf('%d points without a recognisable station code ("sNNN_pMMM") excluded from validation.\n', n_missing);
end
if ~any(has_station)
    error('No point matched the "sNNN_pMMM" station naming pattern: cannot build station means.');
end
station_str = strings(size(name));
station_str(has_station) = string(cellfun(@(c) c{1}, station(has_station), 'UniformOutput', false));

stations = unique(station_str(has_station));
z_ell_station_mean = containers.Map('KeyType', 'char', 'ValueType', 'double');
z_ort_station_mean = containers.Map('KeyType', 'char', 'ValueType', 'double');
n_points_map        = containers.Map('KeyType', 'char', 'ValueType', 'double');
std_within_map       = containers.Map('KeyType', 'char', 'ValueType', 'double');

for k = 1:numel(stations)
    sel = station_str == stations(k);
    z_ell_station_mean(char(stations(k))) = mean(z_ell(sel), 'omitnan');
    z_ort_station_mean(char(stations(k))) = mean(z_ort(sel), 'omitnan');
    n_points_map(char(stations(k)))        = nnz(sel);
    std_within_map(char(stations(k)))       = std(z_ell(sel), 'omitnan');
end
fprintf('%d stations identified from %d points (mean %.1f points/station).\n', ...
    numel(stations), nnz(has_station), mean(cell2mat(values(n_points_map))));

% Broadcast the station mean back to every individual point (X,Y unchanged)
z_ell_mean_pt = nan(size(Xg));
z_ort_mean_pt = nan(size(Xg));
N_points_station   = nan(size(Xg));
std_within_station = nan(size(Xg));
for i = find(has_station)'
    key = char(station_str(i));
    z_ell_mean_pt(i) = z_ell_station_mean(key);
    z_ort_mean_pt(i) = z_ort_station_mean(key);
    N_points_station(i) = n_points_map(key);
    std_within_station(i) = std_within_map(key);
end

% keep only points that belong to a recognised station from here on
keep = has_station;
name = name(keep); Xg = Xg(keep); Yg = Yg(keep);
z_ell_mean_pt = z_ell_mean_pt(keep); z_ort_mean_pt = z_ort_mean_pt(keep);
station_str = station_str(keep);
N_points_station = N_points_station(keep); std_within_station = std_within_station(keep);

%% 3) Bilinear raster sampling at each point's OWN real (X,Y)
[nr, nc] = size(raster_dtm);
[cq, rq] = worldToIntrinsic(R_dtm, Xg, Yg);
in_grid  = isfinite(cq) & isfinite(rq) & cq >= 1 & cq <= nc & rq >= 1 & rq <= nr;
cq(~in_grid) = NaN;  rq(~in_grid) = NaN;

z_orig = samp(raster_dtm,           cq, rq);
z_corr = samp(raster_dtm_corrected, cq, rq);

% membership in the corrected area (nearest pixel) and local DTM roughness (5x5 std)
ri = min(max(round(rq), 1), nr);   ci = min(max(round(cq), 1), nc);
in_corr = false(size(Xg));
in_corr(in_grid) = mask_extended(sub2ind([nr nc], ri(in_grid), ci(in_grid)));
loc_std = nan(size(Xg));
for i = find(in_grid)'
    w = raster_dtm(max(1, ri(i)-2):min(nr, ri(i)+2), max(1, ci(i)-2):min(nc, ci(i)+2));
    loc_std(i) = std(w(:), 'omitnan');
end

%% 4) Automatic vertical datum selection (orthometric vs ellipsoidal)
% Offset computed against the STATION-MEAN elevation, consistently with
% what will be used as the comparison truth throughout.
d_ort = median(z_orig - z_ort_mean_pt, 'omitnan');
d_ell = median(z_orig - z_ell_mean_pt, 'omitnan');
fprintf('Median offset, original DTM - GNSS station mean:  orthometric = %+.2f m | ellipsoidal = %+.2f m\n', d_ort, d_ell);
if abs(d_ell) < abs(d_ort)
    z_gnss_mean = z_ell_mean_pt;  datum = 'ellipsoidal';
else
    z_gnss_mean = z_ort_mean_pt;  datum = 'orthometric';
end
fprintf('Using %s elevation (smaller absolute offset).\n', datum);
if min(abs([d_ort d_ell])) > 5
    warning(['Offset > 5 m with both datums: check the DTM vertical datum/CRS (geoid), ' ...
             'point coordinates, or units.']);
end

%% 5) Residuals and valid-point selection
ok = in_grid & isfinite(z_gnss_mean) & isfinite(z_orig) & isfinite(z_corr);
fprintf('Points in file: %d | valid: %d | off-grid/NaN: %d | inside corrected area: %d\n', ...
    numel(Xg), nnz(ok), nnz(~ok), nnz(ok & in_corr));

e_corr = z_corr - z_gnss_mean;
e_orig = z_orig - z_gnss_mean;

%% 6) Metrics
sets = {ok,           'all'; ...
        ok & in_corr, 'inside corrected area'; ...
        ok & ~in_corr,'outside corrected area'};
M = [];  rn = {};
for s = 1:size(sets, 1)
    sel = sets{s, 1};
    if nnz(sel) == 0, continue; end
    M  = [M; metrics(e_corr(sel)); metrics(e_orig(sel))];               %#ok<AGROW>
    rn = [rn; {['Corrected | ' sets{s,2}]; ['Original | ' sets{s,2}]}]; %#ok<AGROW>
end
Tstats = array2table(M, 'VariableNames', ...
    {'N','ME','Median','SD','RMSE','MAE','NMAD','P95_abs'}, 'RowNames', rn);
disp(Tstats);
fprintf('NB: RMSE^2 = ME^2 + SD^2 (systematic + random error). NMAD = 1.4826*MAD, robust to outliers.\n');
fprintf('N = number of POINTS, each compared against its own station''s mean GNSS elevation.\n');

% Bootstrap confidence intervals (corrected DTM, inside corrected area) and gain vs original
sel = ok & in_corr;
if nnz(sel) >= 10
    ci = boot_ci(e_corr(sel), 2000);
    fprintf('\n95%% bootstrap CI (corrected, inside corrected area):\n');
    fprintf('  ME   [%+.3f, %+.3f] m\n  RMSE [%.3f, %.3f] m\n  MAE  [%.3f, %.3f] m\n  NMAD [%.3f, %.3f] m\n', ...
        ci(1,1), ci(2,1), ci(1,2), ci(2,2), ci(1,3), ci(2,3), ci(1,4), ci(2,4));
    [dR, dCI] = boot_delta_rmse(e_orig(sel), e_corr(sel), 2000);
    fprintf('RMSE reduction (original - corrected): %+.3f m  [95%% CI %+.3f, %+.3f]\n', dR, dCI(1), dCI(2));
    if dCI(1) <= 0
        fprintf('  -> the correction does NOT improve on the original DTM in a statistically distinguishable way.\n');
    end
else
    fprintf('\nFewer than 10 points inside the corrected area: bootstrap skipped.\n');
end

% Note on non-independence: points from the same station share the same
% truth value, so their residuals are correlated by construction. The
% bootstrap above resamples individual points, which slightly overstates
% precision; a stricter alternative is a block-bootstrap by station - see
% the companion point-wise and station-centroid functions for comparison.
n_eff_stations = numel(unique(station_str(ok & in_corr)));
fprintf('(%d points above come from %d distinct stations - residuals are not fully independent.)\n', ...
    nnz(ok & in_corr), n_eff_stations);

%% 7) Per-point table and export
Tpts = table(name, station_str, Xg, Yg, z_gnss_mean, z_orig, z_corr, e_orig, e_corr, ...
    z_orig - z_gnss_mean, z_orig - z_corr, in_corr, loc_std, N_points_station, std_within_station, ok, ...
    'VariableNames', {'PointName','Station','X','Y','z_gnss_station_mean','z_dtm_orig','z_dtm_corr', ...
    'res_orig','res_corr','DoD_observed','DoD_modelled','inside_corrected_area', ...
    'local_std_5x5','N_points_in_station','std_within_station','valid'});

fprintf('\nLargest-residual points (corrected DTM):\n');
tmp = Tpts(ok, :);
[~, ix] = sort(abs(tmp.res_corr), 'descend');
disp(tmp(ix(1:min(5, height(tmp))), {'PointName','Station','X','Y','z_gnss_station_mean','z_dtm_corr','res_corr','inside_corrected_area'}));

if ~isfolder(out_dir), mkdir(out_dir); end
writetable(Tpts,   fullfile(out_dir, 'gnss_validation_mean_records.csv'));
writetable(Tstats, fullfile(out_dir, 'gnss_validation_mean_metrics.csv'), 'WriteRowNames', true);

%% 8) Diagnostic figures
fig = figure('Color', 'w', 'Position', [100 100 1200 900]);
tiledlayout(2, 2, 'TileSpacing', 'compact');
sel = ok;

nexttile; hold on; grid on; axis equal
plot(z_gnss_mean(sel), z_orig(sel), 'o', 'Color', [0.6 0.6 0.6], 'MarkerSize', 4);
plot(z_gnss_mean(sel), z_corr(sel), 'o', 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
lim = [min([z_gnss_mean(sel); z_orig(sel); z_corr(sel)]), max([z_gnss_mean(sel); z_orig(sel); z_corr(sel)])];
plot(lim, lim, 'k--');
xlabel(sprintf('GNSS station-mean elevation (%s) [m]', datum)); ylabel('Model elevation at point [m]');
legend({'Original DTM', 'Corrected DTM', '1:1'}, 'Location', 'best');
title('Elevation: model (per point) vs GNSS station mean');

nexttile; hold on; grid on
histogram(e_orig(sel), 'BinWidth', 0.1, 'FaceAlpha', 0.5);
histogram(e_corr(sel), 'BinWidth', 0.1, 'FaceAlpha', 0.5);
xline(0, 'k-');
xlabel('Residual = model - GNSS station mean [m]'); ylabel('Number of points');
legend({'Original DTM', 'Corrected DTM'}, 'Location', 'best');
title('Residual distribution');

nexttile; hold on; grid on
mm = mean(e_corr(sel));  ss = std(e_corr(sel));
plot((z_corr(sel) + z_gnss_mean(sel)) / 2, e_corr(sel), 'o', 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'none');
yline(mm, 'k-'); yline(mm + 1.96*ss, 'k--'); yline(mm - 1.96*ss, 'k--');
xlabel('Mean of (model, GNSS station mean) [m]'); ylabel('Residual [m]');
title('Bland-Altman (corrected DTM)');

nexttile; hold on; axis equal; grid on
L = min(3, max(abs(e_corr(sel))));
scatter(Xg(sel), Yg(sel), 45, e_corr(sel), 'filled');
n2 = 128;
cmap = [linspace(0.2, 1, n2)', linspace(0.4, 1, n2)', ones(n2, 1); ...
        ones(n2, 1), linspace(1, 0.2, n2)', linspace(1, 0.2, n2)'];
colormap(gca, cmap); caxis([-L L]); cb = colorbar; ylabel(cb, 'Residual [m]');
xlabel('X (UTM 32N) [m]'); ylabel('Y (UTM 32N) [m]');
title('Spatial distribution of residuals (corrected DTM)');

sgtitle('GNSS validation - point-wise, compared against station-mean elevation');

exportgraphics(fig, fullfile(out_dir, 'gnss_validation_mean.png'), 'Resolution', 200);
end

%% ------------------------------------------------------------------------
function z = samp(A, cq, rq)
% Bilinear interpolation; falls back to nearest-neighbour where a bilinear
% neighbour is NaN.
z = interp2(A, cq, rq, 'linear');
bad = isnan(z) & isfinite(cq) & isfinite(rq);
if any(bad)
    z(bad) = interp2(A, cq(bad), rq(bad), 'nearest');
end
end

function m = metrics(r)
% [N ME Median SD RMSE MAE NMAD P95_abs]
r = r(isfinite(r(:)));  r = r(:);
n = numel(r);
if n == 0, m = nan(1, 8); return; end
sa = sort(abs(r));
m = [n, mean(r), median(r), std(r), sqrt(mean(r.^2)), mean(abs(r)), ...
     1.4826 * median(abs(r - median(r))), sa(max(1, ceil(0.95 * n)))];
end

function ci = boot_ci(r, B)
% 95% bootstrap CI of [ME RMSE MAE NMAD]; row 1 = lower bound, row 2 = upper bound.
r = r(isfinite(r(:)));  r = r(:);  n = numel(r);
bs = zeros(B, 4);
for b = 1:B
    rb = r(randi(n, n, 1));
    bs(b, :) = [mean(rb), sqrt(mean(rb.^2)), mean(abs(rb)), 1.4826 * median(abs(rb - median(rb)))];
end
bs = sort(bs);
ci = [bs(max(1, round(0.025*B)), :); bs(min(B, round(0.975*B)), :)];
end

function [d, ci] = boot_delta_rmse(r_orig, r_corr, B)
% Paired bootstrap of the RMSE_original - RMSE_corrected difference
r_orig = r_orig(:);  r_corr = r_corr(:);  n = numel(r_corr);
d = sqrt(mean(r_orig.^2)) - sqrt(mean(r_corr.^2));
dd = zeros(B, 1);
for b = 1:B
    i = randi(n, n, 1);
    dd(b) = sqrt(mean(r_orig(i).^2)) - sqrt(mean(r_corr(i).^2));
end
dd = sort(dd);
ci = [dd(max(1, round(0.025*B))), dd(min(B, round(0.975*B)))];
end
