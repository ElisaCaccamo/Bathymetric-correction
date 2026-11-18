function [Tgroup, Tst] = compare_modules_gnss(T_hm, T_mm, S, all_centerlines, out_dir, opts)
% COMPARE_MODULES_GNSS  Compare the raw DTM and the DTMs corrected by the
% hydraulic (HM) and morphodynamic (MM) modules against the GNSS cross-sections.
%
% Inputs are the per-point tables written by validation_gnss_mean for the two
% modules (gnss_validation_mean_records.csv) and the station table of
% station_reach_class. Residual convention as in validation_gnss_mean:
%   e = DTM sampled at the point - mean GNSS elevation of its cross-section
% Only the valid points inside the corrected area are used.
%
% STATISTICS (one row per group: each reach type, all stations, and the same
% groups without the stations in opts.exclude)
%   RMSE and ME of the raw DTM, HM and MM, with 95% bootstrap intervals
%   obtained by resampling STATIONS (cross-sections), not points: the points
%   of a cross-section share the same reference elevation and are not
%   independent, so resampling points gives intervals that are too narrow.
%   dRMSE = RMSE_MM - RMSE_HM with a paired station bootstrap interval.
%   n_MM_better = stations where the RMSE of MM is lower than that of HM,
%   with the two-sided sign test p-value (paired comparison by station).
%
% FIGURES (out_dir)
%   GNSS_module_comparison.png
%     A  mean residual of every station for raw DTM, HM and MM (stations in
%        network order, background = reach type)
%     B  MM - HM difference at every station [cm], bar colour = module closer to GNSS
%     C  RMSE by group with station-bootstrap intervals, ME as a marker
%   GNSS_station_assignment.png
%     left   map of the stations coloured by reach type (ambiguous ones ringed)
%     right  MM - HM difference vs distance from the upstream node of the segment
%
% opts (optional struct)
%   exclude  cell of stations reported also separately (default {'s006'})
%   B        bootstrap replicates (default 2000)

if nargin < 6, opts = struct(); end
if ~isfield(opts, 'exclude'), opts.exclude = {'s006'}; end
if ~isfield(opts, 'B'), opts.B = 2000; end
if ~isfolder(out_dir), mkdir(out_dir); end
rng(1);

col_raw = [137 135 129] / 255;
col_hm  = [ 42 120 214] / 255;
col_mm  = [235 104  52] / 255;
ink     = [ 43  43  42] / 255;

%% 1) Match the points of the two modules
name_hm = T_hm.PointName;  if ~iscell(name_hm), name_hm = cellstr(name_hm); end
name_mm = T_mm.PointName;  if ~iscell(name_mm), name_mm = cellstr(name_mm); end
st_hm   = T_hm.Station;    if ~iscell(st_hm),   st_hm   = cellstr(st_hm);   end
[tf, loc] = ismember(name_hm, name_mm);
if ~all(tf)
    warning('%d points of the hydraulic table are missing in the morphodynamic table.', nnz(~tf));
end
ok = tf & logical(T_hm.valid) & logical(T_hm.inside_corrected_area);
loc(~tf) = 1;
ok = ok & logical(T_mm.valid(loc)) & logical(T_mm.inside_corrected_area(loc));
d_orig = abs(T_hm.z_dtm_orig(ok) - T_mm.z_dtm_orig(loc(ok)));
if max(d_orig) > 1e-6
    warning('The raw DTM differs between the two tables (max %.3g m): different raw DTM or GNSS datum?', max(d_orig));
end
e0 = T_hm.res_orig;               % raw DTM
e1 = T_hm.res_corr;               % hydraulic module
e2 = nan(size(e1));
e2(tf) = T_mm.res_corr(loc(tf));  % morphodynamic module

%% 2) Per-station sums (fast bootstrap)
stations = S.Station;  if ~iscell(stations), stations = cellstr(stations); end
cls = S.reach_class;   if ~iscell(cls), cls = cellstr(cls); end
ns = numel(stations);
n  = zeros(ns,1);
SS = zeros(ns,3);   % sum of squares  [raw HM MM]
SM = zeros(ns,3);   % sum             [raw HM MM]
for k = 1:ns
    sel = ok & strcmp(st_hm, stations{k});
    E = [e0(sel), e1(sel), e2(sel)];
    n(k) = nnz(sel);
    SS(k,:) = sum(E.^2, 1);
    SM(k,:) = sum(E, 1);
end
rmse_st = sqrt(SS ./ n);
me_st   = SM ./ n;
dz_cm   = 100 * (me_st(:,3) - me_st(:,2));   % mean MM - HM bed difference at the station [cm]
closer  = repmat({'HM'}, ns, 1);
closer(rmse_st(:,3) < rmse_st(:,2)) = {'MM'};
Tst = table(stations, cls, S.segment_code, n, me_st(:,1), me_st(:,2), me_st(:,3), ...
    rmse_st(:,1), rmse_st(:,2), rmse_st(:,3), dz_cm, closer, ...
    'VariableNames', {'Station','reach_class','segment_code','n_points','ME_raw','ME_HM','ME_MM', ...
    'RMSE_raw','RMSE_HM','RMSE_MM','MM_minus_HM_cm','closer'});

%% 3) Groups
class_order = {'single channel','below bifurcation','below confluence','other multi-thread'};
class_order = class_order(ismember(class_order, cls));
excl = ismember(stations, opts.exclude);
gname = {};  gsel = {};
for c = 1:numel(class_order)
    in_c = strcmp(cls, class_order{c});
    gname{end+1} = class_order{c};  gsel{end+1} = in_c; %#ok<AGROW>
    if any(in_c & excl)
        gname{end+1} = [class_order{c} ' without ' strjoin(opts.exclude(:)', ', ')]; %#ok<AGROW>
        gsel{end+1}  = in_c & ~excl; %#ok<AGROW>
    end
end
gname{end+1} = 'all stations';  gsel{end+1} = true(ns,1);
if any(excl)
    gname{end+1} = ['all without ' strjoin(opts.exclude(:)', ', ')];
    gsel{end+1}  = ~excl;
end

ng = numel(gname);
R = nan(ng, 30);
for g = 1:ng
    idx = find(gsel{g});
    m = numel(idx);
    rm = sqrt(sum(SS(idx,:), 1) / sum(n(idx)));
    me = sum(SM(idx,:), 1) / sum(n(idx));
    rb = nan(opts.B, 3);  mb = nan(opts.B, 3);
    for b = 1:opts.B
        j = idx(randi(m, m, 1));
        rb(b,:) = sqrt(sum(SS(j,:), 1) / sum(n(j)));
        mb(b,:) = sum(SM(j,:), 1) / sum(n(j));
    end
    q = @(v) [prctile_simple(v, 2.5), prctile_simple(v, 97.5)];
    dr = rb(:,3) - rb(:,2);
    k_mm = nnz(rmse_st(idx,3) < rmse_st(idx,2));
    k_hm = nnz(rmse_st(idx,3) > rmse_st(idx,2));
    R(g,:) = [m, sum(n(idx)), rm, me, q(rb(:,1)), q(rb(:,2)), q(rb(:,3)), ...
              q(mb(:,1)), q(mb(:,2)), q(mb(:,3)), rm(3) - rm(2), q(dr), ...
              k_mm, k_hm, sign_test(k_mm, k_hm), ...
              mean(dz_cm(idx)), min(dz_cm(idx)), max(dz_cm(idx)), nnz(abs(dz_cm(idx)) > 2)];
end
vn = {'N_stations','N_points','RMSE_raw','RMSE_HM','RMSE_MM','ME_raw','ME_HM','ME_MM', ...
      'RMSE_raw_lo','RMSE_raw_hi','RMSE_HM_lo','RMSE_HM_hi','RMSE_MM_lo','RMSE_MM_hi', ...
      'ME_raw_lo','ME_raw_hi','ME_HM_lo','ME_HM_hi','ME_MM_lo','ME_MM_hi', ...
      'dRMSE_MM_HM','dRMSE_lo','dRMSE_hi','n_MM_better','n_HM_better','p_sign_test', ...
      'dz_mean_cm','dz_min_cm','dz_max_cm','n_abs_dz_over_2cm'};
Tgroup = array2table(R, 'VariableNames', vn);
Tgroup.group = gname(:);

fprintf('\nGNSS comparison raw / hydraulic (HM) / morphodynamic (MM), points inside the corrected area\n');
fprintf('95%% intervals from %d bootstrap replicates of the STATIONS (cross-sections)\n', opts.B);
fprintf('%-38s %3s %4s | %-24s | %-22s | %-22s | %-22s | %s\n', 'group', 'st', 'pts', ...
    'RMSE raw', 'RMSE HM', 'RMSE MM', 'RMSE MM - HM', 'MM better / sign test');
for g = 1:ng
    r = R(g,:);
    fprintf('%-38s %3d %4d | %.2f [%.2f-%.2f] ME %+.2f | %.2f [%.2f-%.2f] %+.2f | %.2f [%.2f-%.2f] %+.2f | %+.3f [%+.3f, %+.3f] | %d/%d  p=%.2g\n', ...
        gname{g}, r(1), r(2), r(3), r(9), r(10), r(6), r(4), r(11), r(12), r(7), ...
        r(5), r(13), r(14), r(8), r(21), r(22), r(23), r(24), r(1), r(26));
end
fprintf('ME after each RMSE interval. MM - HM < 0: the morphodynamic module is closer to GNSS.\n');

writetable(Tst, fullfile(out_dir, 'GNSS_module_comparison_stations.csv'));
writetable(Tgroup, fullfile(out_dir, 'GNSS_module_comparison_groups.csv'));

%% 4) Figure: panels A, B, C
lab = regexprep(stations, '^[sS]', '');
x = (1:ns)';
band = containers.Map({'single channel','below bifurcation','below confluence','other multi-thread'}, ...
    {[0.88 0.87 0.85], [1 1 1], [0.95 0.94 0.92], [0.92 0.91 0.89]});

fig = figure('Color', 'w', 'Position', [60 60 1400 1050]);

% --- A ---
axA = axes('Position', [0.07 0.60 0.91 0.29]);  hold(axA, 'on');
yl = [min([me_st(:); -0.2]) - 0.1, max(me_st(:)) + 0.15];
hb = draw_bands(axA, cls, yl, band, class_order);
plot(axA, [0.4 ns + 0.6], [0 0], '-', 'Color', ink, 'LineWidth', 1);
h0 = plot(axA, x, me_st(:,1), 'o', 'MarkerSize', 7, 'LineWidth', 1.4, 'Color', col_raw, 'MarkerFaceColor', 'w');
h1 = plot(axA, x - 0.15, me_st(:,2), 'o', 'MarkerSize', 6.5, 'MarkerFaceColor', col_hm, 'MarkerEdgeColor', 'w');
h2 = plot(axA, x + 0.15, me_st(:,3), 'o', 'MarkerSize', 6.5, 'MarkerFaceColor', col_mm, 'MarkerEdgeColor', 'w');
set(axA, 'XLim', [0.4 ns + 0.6], 'YLim', yl, 'XTick', x, 'XTickLabel', lab, 'FontSize', 8, 'Box', 'off', 'Layer', 'top');
grid(axA, 'on');  set(axA, 'XGrid', 'off');
ylabel(axA, {'DTM - GNSS station mean [m]', '> 0: DTM too high'});
title(axA, 'A. Mean residual at each GNSS cross-section (stations in network order, upstream to downstream)', ...
    'FontWeight', 'normal');
legend(axA, [h0 h1 h2 hb], [{'raw LiDAR DTM', 'hydraulic module', 'morphodynamic module'}, class_order], ...
    'Orientation', 'horizontal', 'Location', 'northoutside', 'Box', 'off', 'FontSize', 8, 'NumColumns', 4);

% --- B ---
axB = axes('Position', [0.07 0.37 0.91 0.14]);  hold(axB, 'on');
ylb = [min([dz_cm; -1]) - 1, max([dz_cm; 1]) + 1];
draw_bands(axB, cls, ylb, band, class_order);
is_mm = strcmp(closer, 'MM');
for k = 1:ns
    if is_mm(k), fc = col_mm; else, fc = col_hm; end
    rect(axB, x(k), dz_cm(k), 0.7, fc, 'none', 0.5);
end
plot(axB, [0.4 ns + 0.6], [0 0], '-', 'Color', ink, 'LineWidth', 0.8);
set(axB, 'XLim', [0.4 ns + 0.6], 'YLim', ylb, 'XTick', x, 'XTickLabel', lab, 'FontSize', 8, 'Box', 'off', 'Layer', 'top');
grid(axB, 'on');  set(axB, 'XGrid', 'off');
ylabel(axB, {'MM - HM [cm]'});
title(axB, 'B. Difference between the two corrected DTMs (bar colour = module closer to GNSS: blue hydraulic, orange morphodynamic)', ...
    'FontWeight', 'normal');

% --- C ---
axC = axes('Position', [0.07 0.06 0.91 0.19]);  hold(axC, 'on');
w = 0.26;
cols = {col_raw, col_hm, col_mm};
hC = gobjects(1, 4);
for g = 1:ng
    for j = 1:3
        xx = g + (j - 2) * w;
        rm = R(g, 2 + j);
        lo = R(g, 7 + 2*j);  hi = R(g, 8 + 2*j);
        if j == 1
            hC(j) = rect(axC, xx, rm, w * 0.88, [1 1 1], cols{j}, 1.4);
        else
            hC(j) = rect(axC, xx, rm, w * 0.88, cols{j}, 'none', 0.5);
        end
        plot(axC, [xx xx], [lo hi], '-', 'Color', ink, 'LineWidth', 1);
        plot(axC, xx + [-0.05 0.05], [lo lo], '-', 'Color', ink, 'LineWidth', 1);
        plot(axC, xx + [-0.05 0.05], [hi hi], '-', 'Color', ink, 'LineWidth', 1);
        hC(4) = plot(axC, xx, R(g, 5 + j), 'd', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ink);
        text(axC, xx, hi + 0.02, sprintf('%.2f', rm), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 7, 'Color', ink);
    end
end
plot(axC, [0.5 ng + 0.5], [0 0], '-', 'Color', ink, 'LineWidth', 0.8);
glab = cellfun(@(s, k) sprintf('%s (%d)', s, k), gname(:), num2cell(R(:,1)), 'UniformOutput', false);
set(axC, 'XLim', [0.5 ng + 0.5], 'XTick', 1:ng, 'XTickLabel', glab, 'FontSize', 8, 'Box', 'off');
yc = get(axC, 'YLim');  set(axC, 'YLim', [min(yc(1), 0) max(max(R(:, [10 12 14]))) + 0.12]);
grid(axC, 'on');  set(axC, 'XGrid', 'off');
ylabel(axC, 'RMSE, ME [m]');
title(axC, 'C. RMSE by reach type (number of stations) with 95% station-bootstrap intervals; diamond = mean error ME', ...
    'FontWeight', 'normal');
legend(axC, hC, {'raw LiDAR DTM', 'hydraulic module', 'morphodynamic module', 'ME (mean error)'}, ...
    'Orientation', 'horizontal', 'Location', 'northoutside', 'Box', 'off', 'FontSize', 8);

exportgraphics(fig, fullfile(out_dir, 'GNSS_module_comparison.png'), 'Resolution', 200);

%% 5) Figure: station assignment check
mk = containers.Map({'single channel','below bifurcation','below confluence','other multi-thread'}, ...
    {'s', '^', 'd', 'o'});
fig2 = figure('Color', 'w', 'Position', [80 80 1400 650]);
ax1 = axes('Position', [0.05 0.1 0.5 0.82]);  hold(ax1, 'on');
plot(ax1, all_centerlines(:,2), all_centerlines(:,3), '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 4);
codes = unique(all_centerlines(:,6));
for c = reshape(codes, 1, [])
    r = find(all_centerlines(:,6) == c);
    rm_ = r(ceil(end/2));
    text(ax1, all_centerlines(rm_,2), all_centerlines(rm_,3), sprintf(' %d', c), 'FontSize', 7, 'Color', [0.45 0.45 0.45]);
end
hm = gobjects(1, numel(class_order));
for c = 1:numel(class_order)
    in_c = strcmp(cls, class_order{c});
    hm(c) = plot(ax1, S.X(in_c), S.Y(in_c), mk(class_order{c}), 'MarkerSize', 7, ...
        'MarkerFaceColor', ink, 'MarkerEdgeColor', 'w');
end
amb = logical(S.ambiguous);
ha = plot(ax1, S.X(amb), S.Y(amb), 'o', 'MarkerSize', 14, 'Color', [0.8 0.1 0.1], 'LineWidth', 1.3);
text(ax1, S.X, S.Y, strcat({'  '}, lab), 'FontSize', 7, 'Color', ink);
axis(ax1, 'equal');  grid(ax1, 'on');
xlabel(ax1, 'X (UTM 32N) [m]');  ylabel(ax1, 'Y (UTM 32N) [m]');
title(ax1, 'Stations by reach type (grey numbers = segment codes)', 'FontWeight', 'normal');
if any(amb)
    legend(ax1, [hm ha], [class_order, {'ambiguous (points on two segments)'}], 'Location', 'best', 'Box', 'off', 'FontSize', 8);
else
    legend(ax1, hm, class_order, 'Location', 'best', 'Box', 'off', 'FontSize', 8);
end

ax2 = axes('Position', [0.62 0.1 0.36 0.82]);  hold(ax2, 'on');
for c = 1:numel(class_order)
    in_c = strcmp(cls, class_order{c});
    plot(ax2, S.s_from_node(in_c), dz_cm(in_c), mk(class_order{c}), 'MarkerSize', 7, ...
        'MarkerFaceColor', ink, 'MarkerEdgeColor', 'w');
end
plot(ax2, [0 max(S.s_from_node) * 1.05 + 1], [0 0], '-', 'Color', ink);
text(ax2, S.s_from_node, dz_cm, strcat({' '}, lab), 'FontSize', 7, 'Color', ink);
grid(ax2, 'on');
xlabel(ax2, 'Distance from the upstream node of the segment [m]');
ylabel(ax2, 'MM - HM at the station [cm]');
title(ax2, 'Effect of the morphodynamic module vs distance from the node', 'FontWeight', 'normal');

exportgraphics(fig2, fullfile(out_dir, 'GNSS_station_assignment.png'), 'Resolution', 200);
end

%% ------------------------------------------------------------------------
function h = draw_bands(ax, cls, yl, band, class_order)
% Background rectangles over consecutive stations of the same reach type.
% Returns one legend patch per reach type (invisible, with a grey outline so
% that the white band is visible in the legend).
i = 1;  ns = numel(cls);
while i <= ns
    j = i;
    while j < ns && strcmp(cls{j+1}, cls{i}), j = j + 1; end
    patch(ax, [i-0.5 j+0.5 j+0.5 i-0.5], [yl(1) yl(1) yl(2) yl(2)], band(cls{i}), 'EdgeColor', 'none');
    i = j + 1;
end
h = gobjects(1, numel(class_order));
for c = 1:numel(class_order)
    h(c) = patch(ax, NaN(1,4), NaN(1,4), band(class_order{c}), 'EdgeColor', [0.75 0.75 0.75]);
end
end

function h = rect(ax, x, y, w, face, edge, lw)
% Bar drawn as a patch from 0 to y, centred on x, width w (same look in every MATLAB version).
h = patch(ax, x + w/2 * [-1 1 1 -1], [0 0 y y], face, 'EdgeColor', edge, 'LineWidth', lw);
end

function v = prctile_simple(x, p)
% Percentile without the Statistics Toolbox (linear interpolation).
x = sort(x(isfinite(x)));
if isempty(x), v = NaN; return; end
pos = 1 + (numel(x) - 1) * p / 100;
lo = floor(pos);  hi = ceil(pos);
v = x(lo) + (pos - lo) * (x(hi) - x(lo));
end

function p = sign_test(k1, k2)
% Two-sided exact sign test (ties already excluded).
nn = k1 + k2;
if nn == 0, p = NaN; return; end
k = min(k1, k2);
i = 0:k;
logc = gammaln(nn + 1) - gammaln(i + 1) - gammaln(nn - i + 1) - nn * log(2);
p = min(1, 2 * sum(exp(logc)));
end
