function [P, forced] = read_field_points(gnss_path, adcp_path, station_split, adcp_min_spacing)
% READ_FIELD_POINTS  Field points of the validation in one table: GNSS
% cross-sections surveyed by wading and ADCP transects. The two surveys are
% treated in the same way: every cross-section (GNSS) or transect (ADCP) is a
% "station", compared with its mean bed elevation (validation_gnss_mean).
%
% GNSS file  csv "NomePunto, Y, X, Quota ellissoidica, Quota ortometrica",
%            point names "sNNN_pMMM" -> station "sNNN", point number MMM
% ADCP file  csv "N_UTM32N, E_UTM32N, Z_BOTTOM" (+ a 4th field that is filled
%            only on the first row of each transect, e.g. "20250307151356_-1sec").
%            Each transect becomes a station "aNN" (a01, a02, ... in file
%            order), points "aNN_pMMM". Z_BOTTOM is an ellipsoidal elevation:
%            checked on 19/11/2026 against the ellipsoidal elevation of the
%            nearest GNSS points (4-9 cm apart; the orthometric one is ~48 m
%            lower), i.e. the same reference as the DTM. z_ort is NaN.
%
% station_split  (optional) cell n x 4 {station, first point, last point,
%            segment code}: GNSS cross-sections that cross more than one
%            segment are split into sub-stations, one per segment, named
%            "sNNN-tCC" (CC = segment code; hyphen, not underscore, so that the
%            name is not rendered as a TeX subscript in the figures). Each sub-station is compared with
%            the mean of its own points, and its segment is the one given here
%            (manual assignment, see station_reach_class).
% adcp_min_spacing  (optional, default 0) minimum distance [m] between two
%            kept ADCP points along a transect; 0 keeps all the points.
%
% OUTPUTS
%   P       table: PointName, Station, X, Y, z_ell, z_ort, source ('GNSS' or 'ADCP')
%   forced  table: Station, segment_code of the sub-stations of station_split

if nargin < 2, adcp_path = ''; end
if nargin < 3 || isempty(station_split), station_split = cell(0, 4); end
if nargin < 4 || isempty(adcp_min_spacing), adcp_min_spacing = 0; end

%% GNSS (same parsing as before in validation_gnss_mean)
opts = detectImportOptions(gnss_path, 'VariableNamingRule', 'preserve');
T = readtable(gnss_path, opts);
if width(T) < 5
    error('GNSS file: expected >= 5 columns, found %d. Check the delimiter.', width(T));
end
V = nan(height(T), 4);
for k = 1:4
    c = T{:, k+1};
    if isnumeric(c)
        V(:, k) = double(c);
    else
        V(:, k) = str2double(strrep(string(c), ',', '.'));   % comma decimal mark
    end
end
name = string(T{:, 1});
Yg = V(:, 1);  Xg = V(:, 2);   % NB: Y precedes X in the source file
if median(Xg, 'omitnan') > median(Yg, 'omitnan')   % UTM 32N: Easting ~4e5 < Northing ~5e6
    warning('GNSS file: X and Y look swapped relative to UTM32N expectations: swapping them back.');
    [Xg, Yg] = deal(Yg, Xg);
end
tok = regexp(name, '^(s\d+)_p(\d+)$', 'tokens', 'once', 'ignorecase');
has = ~cellfun(@isempty, tok);
if nnz(~has) > 0
    fprintf('%d GNSS points without a recognisable station code ("sNNN_pMMM") excluded from validation.\n', nnz(~has));
end
station = strings(size(name));
pnum = nan(size(name));
station(has) = string(cellfun(@(c) lower(c{1}), tok(has), 'UniformOutput', false));
pnum(has) = cellfun(@(c) str2double(c{2}), tok(has));

%% Split of the GNSS cross-sections that cross more than one segment
forced = table(strings(0,1), zeros(0,1), 'VariableNames', {'Station','segment_code'});
split_st = unique(string(station_split(:,1)));
new_station = station;
for r = 1:size(station_split, 1)
    st = lower(string(station_split{r,1}));
    p1 = station_split{r,2};  p2 = station_split{r,3};  code = station_split{r,4};
    sub = sprintf('%s-t%02d', st, code);
    sel = station == st & pnum >= p1 & pnum <= p2;
    if ~any(sel)
        warning('Station split: no point of %s between p%03d and p%03d.', st, p1, p2);
    end
    if any(new_station(sel) ~= station(sel))
        warning('Station split: points of %s assigned to more than one sub-station.', st);
    end
    new_station(sel) = sub;
    forced = [forced; {string(sub), code}]; %#ok<AGROW>
end
for st = reshape(split_st, 1, [])
    left = new_station == lower(st);
    if any(left)
        warning('Station split: %d points of %s are not in any sub-station and are excluded.', nnz(left), st);
        new_station(left) = "";
    end
end
station = new_station;

Pg = table(name(has), station(has), Xg(has), Yg(has), V(has,3), V(has,4), ...
    repmat("GNSS", nnz(has), 1), ...
    'VariableNames', {'PointName','Station','X','Y','z_ell','z_ort','source'});
Pg = Pg(strlength(Pg.Station) > 0, :);

%% ADCP
Pa = Pg([], :);
if ~isempty(adcp_path)
    txt = fileread(adcp_path);
    lines = regexp(txt, '\r\n|\r|\n', 'split');
    lines = lines(~cellfun(@(s) isempty(strtrim(s)), lines));
    N = []; E = []; Z = []; tr = []; k_tr = 0;
    for i = 1:numel(lines)
        f = strsplit(lines{i}, ',', 'CollapseDelimiters', false);
        if numel(f) < 3, continue; end
        v = str2double(f(1:3));
        if any(isnan(v)), continue; end          % header or malformed row
        if (numel(f) >= 4 && ~isempty(strtrim(f{4}))) || k_tr == 0
            k_tr = k_tr + 1;                     % first row of a new transect
        end
        N(end+1,1) = v(1); E(end+1,1) = v(2); Z(end+1,1) = v(3); tr(end+1,1) = k_tr; %#ok<AGROW>
    end
    keep = true(size(N));
    if adcp_min_spacing > 0
        for t = unique(tr)'
            idx = find(tr == t);
            last = idx(1);
            for j = idx(2:end)'
                if hypot(E(j) - E(last), N(j) - N(last)) < adcp_min_spacing
                    keep(j) = false;
                else
                    last = j;
                end
            end
        end
    end
    pname = strings(size(N));  st_a = strings(size(N));
    for t = unique(tr)'
        idx = find(tr == t);
        st_a(idx) = sprintf('a%02d', t);
        pname(idx) = string(sprintf('a%02d_p', t)) + compose('%03d', (1:numel(idx))');
    end
    Pa = table(pname(keep), st_a(keep), E(keep), N(keep), Z(keep), nan(nnz(keep),1), ...
        repmat("ADCP", nnz(keep), 1), ...
        'VariableNames', {'PointName','Station','X','Y','z_ell','z_ort','source'});
    fprintf('ADCP: %d transects, %d points (%d kept, minimum spacing %.2f m).\n', ...
        max([tr; 0]), numel(N), nnz(keep), adcp_min_spacing);
end

P = [Pg; Pa];
fprintf('Field points: %d GNSS in %d stations, %d ADCP in %d transects.\n', ...
    height(Pg), numel(unique(Pg.Station)), height(Pa), numel(unique(Pa.Station)));
end
