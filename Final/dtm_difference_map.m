function [dZ, Tstats] = dtm_difference_map(dtm_hm, dtm_mm, dtm_raw, R, all_centerlines, S, out_dir, crs_code)
% DTM_DIFFERENCE_MAP  Difference between the DTMs corrected by the
% morphodynamic (MM) and hydraulic (HM) modules, dZ = DTM_MM - DTM_HM.
%
%   dZ < 0  the morphodynamic bed is LOWER (deeper) than the hydraulic one
%   dZ > 0  the morphodynamic bed is HIGHER than the hydraulic one
% dZ is computed only where at least one of the two modules changed the raw
% DTM; elsewhere it is NaN. Both corrected DTMs are capped at the raw DTM,
% so where the raw DTM is the cap of both, dZ = 0.
%
% It shows where the morphodynamic module really acts. The morphodynamic
% module changes the split of the discharge at the bifurcations, so it changes
% the bed in every segment downstream of the first bifurcation, not only in
% the first reach of the branches. In a single channel it should give dZ = 0
% (same discharge, width and slope): any difference there points to a
% different discharge in the two modules (see check_discharge_conservation).
%
% INPUTS
%   dtm_hm, dtm_mm, dtm_raw  corrected (HM, MM) and raw DTM on the same grid (NaN = nodata)
%   R                        spatial reference of the grid
%   all_centerlines          [ID X Y z_DTM W segment_code] (segment labels)
%   S                        station table of station_reach_class (station labels), or []
%   out_dir                  output folder
%   crs_code                 EPSG code of the output GeoTIFF
%
% OUTPUTS
%   dZ      difference raster (also written as DTM_difference_MM_minus_HM.tif)
%   Tstats  summary of dZ over the corrected pixels

if ~isfolder(out_dir), mkdir(out_dir); end
changed = isfinite(dtm_raw) & isfinite(dtm_hm) & isfinite(dtm_mm) & ...
          (abs(dtm_hm - dtm_raw) > 1e-6 | abs(dtm_mm - dtm_raw) > 1e-6);
dZ = dtm_mm - dtm_hm;
dZ(~changed) = NaN;

v = dZ(changed);
n = numel(v);
Tstats = table(n, mean(v), median(v), min(v), max(v), ...
    nnz(abs(v) <= 0.005) / n, nnz(v < -0.005) / n, nnz(v > 0.005) / n, ...
    prctile_simple(abs(v), 95), ...
    'VariableNames', {'N_pixels','dZ_mean','dZ_median','dZ_min','dZ_max', ...
    'share_abs_le_5mm','share_MM_deeper','share_MM_higher','P95_abs'});
fprintf('\nDTM_MM - DTM_HM over %d corrected pixels: mean %+.3f m, median %+.3f m, range [%+.3f, %+.3f] m\n', ...
    n, Tstats.dZ_mean, Tstats.dZ_median, Tstats.dZ_min, Tstats.dZ_max);
fprintf('  |dZ| <= 5 mm: %.1f%% | MM deeper: %.1f%% | MM higher: %.1f%% | P95 |dZ| = %.3f m\n', ...
    100 * Tstats.share_abs_le_5mm, 100 * Tstats.share_MM_deeper, 100 * Tstats.share_MM_higher, Tstats.P95_abs);
writetable(Tstats, fullfile(out_dir, 'DTM_difference_MM_minus_HM_stats.csv'));
geotiffwrite(fullfile(out_dir, 'DTM_difference_MM_minus_HM.tif'), dZ, R, 'CoordRefSysCode', crs_code);

%% Map
% Diverging colour scale centred on zero: blue = MM deeper, red = MM higher,
% grey = no difference. Symmetric limits at the 98th percentile of |dZ|.
L = max(prctile_simple(abs(v), 98), 0.01);
cneg = [0.16 0.37 0.70];  cmid = [0.92 0.92 0.91];  cpos = [0.75 0.22 0.10];
m = 128;  t = linspace(0, 1, m)';
cmap = [cneg + t .* (cmid - cneg); cmid + t .* (cpos - cmid)];

[nr, nc] = size(dZ);
[x_img, ~] = intrinsicToWorld(R, [1 nc], [1 1]);
[~, y_img] = intrinsicToWorld(R, [1 1], [1 nr]);

fig = figure('Color', 'w', 'Position', [60 60 1300 900]);
ax = axes('Position', [0.06 0.07 0.82 0.86]);
imagesc(ax, x_img, y_img, dZ, 'AlphaData', double(isfinite(dZ)));
axis(ax, 'xy');  axis(ax, 'image');  hold(ax, 'on');
set(ax, 'Color', 'w');
colormap(ax, cmap);  caxis(ax, [-L L]);
cb = colorbar(ax);
ylabel(cb, 'DTM_{MM} - DTM_{HM} [m]   (< 0: morphodynamic bed deeper)');

codes = unique(all_centerlines(:,6));
for c = reshape(codes, 1, [])
    r = find(all_centerlines(:,6) == c);
    rm = r(ceil(end/2));
    text(ax, all_centerlines(rm,2), all_centerlines(rm,3), sprintf('%d', c), 'FontSize', 8, ...
        'Color', [0.25 0.25 0.25], 'FontWeight', 'bold');
end
if ~isempty(S)
    plot(ax, S.X, S.Y, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w');
    text(ax, S.X, S.Y, strcat({'  '}, S.Station), 'FontSize', 7, 'Color', 'k');
end
xlabel(ax, 'X (UTM 32N) [m]');  ylabel(ax, 'Y (UTM 32N) [m]');
title(ax, {'Difference between the DTMs corrected by the morphodynamic and hydraulic modules', ...
    'bold numbers = segment codes, black dots = in situ stations (GNSS, ADCP)'}, 'FontWeight', 'normal');
exportgraphics(fig, fullfile(out_dir, 'DTM_difference_MM_minus_HM.png'), 'Resolution', 200);
end

%% ------------------------------------------------------------------------
function v = prctile_simple(x, p)
% Percentile without the Statistics Toolbox (linear interpolation).
x = sort(x(isfinite(x)));
if isempty(x), v = NaN; return; end
pos = 1 + (numel(x) - 1) * p / 100;
lo = floor(pos);  hi = ceil(pos);
v = x(lo) + (pos - lo) * (x(hi) - x(lo));
end
