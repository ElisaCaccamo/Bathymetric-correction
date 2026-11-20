function Tsec = gnss_cross_section_profiles(gnss_path, T_hm, T_mm, stations, S, all_centerlines, ws_hm, ws_mm, dtm_raw, R, out_dir)
% GNSS_CROSS_SECTION_PROFILES  Transverse profiles of selected GNSS stations:
% GNSS elevation of every point, raw DTM, DTM corrected by the hydraulic (HM)
% and morphodynamic (MM) modules, plus the model data of the nearest sections.
%
% The validation tables keep only the station-mean GNSS elevation, so the
% elevation of every single point is read again from the GNSS file. The
% vertical datum (ellipsoidal/orthometric) is the one whose station mean
% matches z_gnss_station_mean of the validation table.
%
% For every station the figure shows, along the cross-section (distance
% from one end of the measured section, along its main direction):
%   black dots   GNSS elevation of each point;  black dashed line = station mean
%   grey line    raw LiDAR DTM at the points
%   blue / orange   corrected DTM of HM / MM at the points
%   dotted grey  z_DTM of the reference section (water surface used by the modules)
%   thin dashed blue / orange   flat bed of the reference section in HM / MM
% Hollow markers: points outside the corrected area.
%
% The table (printed and written as csv) lists, for the reference section and
% the 3 sections upstream and downstream of the same segment:
%   z_DTM (water surface of the centerline), raw DTM sampled at the centerline
%   point, width W, slope S, and Q, h, z of both modules.
% It is meant to check whether the depth imposed by the modules is plausible
% (e.g. a depth of ~2 m where the GNSS shows ~0.2 m of water) and whether
% z_DTM oscillates along the segment.
%
% INPUTS
%   gnss_path     GNSS csv ("NomePunto, Y, X, Quota ellissoidica, Quota ortometrica")
%   T_hm, T_mm    per-point tables of validation_gnss_mean (HM, MM)
%   stations      cell of station names, e.g. {'s006','s019','s020'}
%   S             station table of station_reach_class
%   all_centerlines  [ID X Y z_DTM W segment_code]
%   ws_hm         struct with fields Q_map, h_rect_all, z_rect_all, widths_all, slope_sect
%   ws_mm         struct with fields Q_map, h_rect_morpho, z_rect_morpho
%   dtm_raw, R    raw DTM and its spatial reference
%   out_dir       output folder

if ~isfolder(out_dir), mkdir(out_dir); end
col_raw = [137 135 129] / 255;
col_hm  = [ 42 120 214] / 255;
col_mm  = [235 104  52] / 255;

%% GNSS file (same parsing as validation_gnss_mean)
opts = detectImportOptions(gnss_path, 'VariableNamingRule', 'preserve');
T = readtable(gnss_path, opts);
V = nan(height(T), 4);
for k = 1:4
    c = T{:, k+1};
    if isnumeric(c), V(:,k) = double(c); else, V(:,k) = str2double(strrep(cellstr(c), ',', '.')); end
end
gname = T{:,1};  if ~iscell(gname), gname = cellstr(gname); end
z_ell = V(:,3);  z_ort = V(:,4);

name_hm = T_hm.PointName;  if ~iscell(name_hm), name_hm = cellstr(name_hm); end
name_mm = T_mm.PointName;  if ~iscell(name_mm), name_mm = cellstr(name_mm); end
st_hm = T_hm.Station;      if ~iscell(st_hm), st_hm = cellstr(st_hm); end
S_st = S.Station;          if ~iscell(S_st), S_st = cellstr(S_st); end

ns = numel(stations);
fig = figure('Color', 'w', 'Position', [60 60 1100 330 * ns]);
rows_out = {};
for k = 1:ns
    st = stations{k};
    sel = find(strcmp(st_hm, st) & logical(T_hm.valid));
    if isempty(sel)
        warning('Station %s not found in the validation table.', st);
        continue
    end
    [~, ig] = ismember(lower(name_hm(sel)), lower(gname));
    [~, im] = ismember(name_hm(sel), name_mm);
    zm = T_hm.z_gnss_station_mean(sel(1));
    if ismember('z_point', T_hm.Properties.VariableNames)
        zg = T_hm.z_point(sel);   % own elevation of each point (GNSS or ADCP), datum of the validation
        datum = 'datum of the validation';
    else
        good = ig > 0;
        if abs(mean(z_ell(ig(good))) - zm) <= abs(mean(z_ort(ig(good))) - zm)
            zg_all = z_ell;  datum = 'ellipsoidal';
        else
            zg_all = z_ort;  datum = 'orthometric';
        end
        zg = nan(numel(sel), 1);  zg(good) = zg_all(ig(good));
    end

    % distance along the main direction of the measured section
    P = [T_hm.X(sel), T_hm.Y(sel)];
    P0 = P - mean(P, 1);
    [~, ~, Vv] = svd(P0, 'econ');
    d = P0 * Vv(:,1);
    d = d - min(d);
    [d, o] = sort(d);
    sel = sel(o);  zg = zg(o);  im = im(o);
    inside = logical(T_hm.inside_corrected_area(sel));
    z_mm_pts = nan(numel(sel), 1);
    z_mm_pts(im > 0) = T_mm.z_dtm_corr(im(im > 0));

    % reference section of the station
    j = find(strcmp(S_st, st), 1);
    if isempty(j)
        row = [];  ttl = st;
    else
        row = S.section_ID(j);
        rc = S.reach_class;  if ~iscell(rc), rc = cellstr(rc); end
        ttl = sprintf('%s - segment %d, section %d (%s)', st, S.segment_code(j), row, rc{j});
    end

    ax = subplot(ns, 1, k);  hold(ax, 'on');
    plot(ax, d, T_hm.z_dtm_orig(sel), '-', 'Color', col_raw, 'LineWidth', 1.8);
    plot(ax, d, T_hm.z_dtm_corr(sel), '-', 'Color', col_hm, 'LineWidth', 1.8);
    plot(ax, d, z_mm_pts, '-', 'Color', col_mm, 'LineWidth', 1.8);
    plot(ax, d, zg, 'k-', 'LineWidth', 0.8);
    plot(ax, d(inside), zg(inside), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    plot(ax, d(~inside), zg(~inside), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
    plot(ax, [d(1) d(end)], [zm zm], 'k--', 'LineWidth', 1);
    if ~isempty(row)
        plot(ax, [d(1) d(end)], all_centerlines(row,4) * [1 1], ':', 'Color', col_raw, 'LineWidth', 1.5);
        plot(ax, [d(1) d(end)], ws_hm.z_rect_all(row) * [1 1], '--', 'Color', col_hm, 'LineWidth', 0.8);
        plot(ax, [d(1) d(end)], ws_mm.z_rect_morpho(row) * [1 1], '--', 'Color', col_mm, 'LineWidth', 0.8);
    end
    grid(ax, 'on');  box(ax, 'off');
    xlabel(ax, 'Distance along the measured cross-section [m]');
    ylabel(ax, sprintf('Elevation (%s) [m]', datum));
    title(ax, ttl, 'FontWeight', 'normal', 'Interpreter', 'none');
    if k == 1
        legend(ax, {'raw LiDAR DTM', 'DTM hydraulic module', 'DTM morphodynamic module', ...
            'in situ points', 'in situ (inside corrected area)', 'in situ (outside)', 'in situ station mean', ...
            'z_DTM of the section (water surface)', 'flat bed HM', 'flat bed MM'}, ...
            'Location', 'eastoutside', 'Box', 'off', 'FontSize', 8, 'Interpreter', 'none');
    end
    fprintf('\n%s: %d points, GNSS %s, station mean %.3f m, raw DTM - GNSS mean = %+.2f m, HM %+.2f m, MM %+.2f m\n', ...
        st, numel(sel), datum, zm, mean(T_hm.z_dtm_orig(sel) - zm), mean(T_hm.z_dtm_corr(sel) - zm), ...
        mean(z_mm_pts - zm, 'omitnan'));
    fprintf('   GNSS: min %.3f, max %.3f m (relief %.2f m) | raw DTM: min %.3f, max %.3f m (relief %.2f m)\n', ...
        min(zg), max(zg), max(zg) - min(zg), min(T_hm.z_dtm_orig(sel)), max(T_hm.z_dtm_orig(sel)), ...
        max(T_hm.z_dtm_orig(sel)) - min(T_hm.z_dtm_orig(sel)));

    % model data of the neighbouring sections of the same segment
    if ~isempty(row)
        seg = all_centerlines(row, 6);
        rr = (row - 3):(row + 3);
        rr = rr(rr >= 1 & rr <= size(all_centerlines, 1));
        rr = rr(all_centerlines(rr, 6) == seg);
        [ci, ri] = worldToIntrinsic(R, all_centerlines(rr,2), all_centerlines(rr,3));
        z_raw_cl = interp2(dtm_raw, ci, ri, 'nearest');
        for q = rr(:)'
            rows_out(end+1, :) = {st, q, seg, all_centerlines(q,4), z_raw_cl(rr == q), ...
                ws_hm.widths_all(q), ws_hm.slope_sect(q), ...
                ws_hm.Q_map(q), ws_hm.h_rect_all(q), ws_hm.z_rect_all(q), ...
                ws_mm.Q_map(q), ws_mm.h_rect_morpho(q), ws_mm.z_rect_morpho(q)}; %#ok<AGROW>
        end
    end
end
exportgraphics(fig, fullfile(out_dir, 'GNSS_cross_section_profiles.png'), 'Resolution', 200);

vn = {'Station','section_ID','segment_code','z_DTM_centerline','z_raw_DTM_at_centerline', ...
      'W','S','Q_HM','h_HM','z_bed_HM','Q_MM','h_MM','z_bed_MM'};
if isempty(rows_out)
    Tsec = cell2table(cell(0, numel(vn)), 'VariableNames', vn);
else
    Tsec = cell2table(rows_out, 'VariableNames', vn);
end
fprintf('\nModel data of the sections around the selected stations:\n');
disp(Tsec);
writetable(Tsec, fullfile(out_dir, 'GNSS_cross_section_profiles_sections.csv'));
end
