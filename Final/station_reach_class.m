function S = station_reach_class(Tpts, all_centerlines, sections_ids, predecessor, bif_info, conf_info, network_root, topo_order, no_flow_segments)
% STATION_REACH_CLASS  Assign every GNSS station to a river segment and to a
% reach type of the network.
%
% A station is one measured GNSS cross-section: all the points named
% "sNNN_pMMM" with the same "sNNN". Only the valid points inside the
% corrected area are used.
%
% Steps
%   1) every point is assigned to the nearest centerline section; the station
%      takes the segment of the majority of its points. share_main is the
%      fraction of points on that segment: share_main < 1 means that the
%      cross-section touches two segments (ambiguous assignment);
%   2) the reference section of the station is the section of that segment
%      nearest to the station centroid (dist_axis); dist_other is the distance
%      from the centroid to the nearest section of any OTHER segment;
%   3) s_from_node is the distance along the centerline from the first
%      section of the segment (a node: bifurcation, confluence or network
%      head) to the reference section;
%   4) the flow fraction of each segment is a topological count of threads:
%      1 at the network head, halved in each branch of a bifurcation, summed
%      at a confluence (it is not a discharge). A segment without discharge
%      at the flight (no_flow_segments) has fraction 0, and the other branch
%      of its bifurcation keeps the whole fraction of the mother: that node is
%      not a real split at the flight;
%   5) reach type of the segment:
%        'single channel'      flow fraction = 1
%        'below bifurcation'   segment starting at a bifurcation
%        'below confluence'    multi-thread segment starting at a confluence
%        'other multi-thread'  any other multi-thread segment
%        'no flow'             segment in no_flow_segments
%
% INPUTS
%   Tpts             per-point table of validation_gnss_mean
%                    (gnss_validation_mean_records.csv)
%   all_centerlines  [ID X Y z_DTM W segment_code], ID = row index
%   sections_ids     [first_ID last_ID] of each segment (row = segment index)
%   predecessor      cell, upstream segment index(es) of each segment
%   bif_info         cell {mother_ID, [first_ID_branch1 first_ID_branch2]}
%   conf_info        cell {confluence_ID (first section downstream), [last_IDs upstream]}
%   network_root     segment index of the network head
%   topo_order       segment indices in topological order (upstream first)
%   no_flow_segments (optional) codes of the segments without discharge at the flight
%
% OUTPUT
%   S  table, one row per station, sorted along the network (topological
%      order of the segments, then section ID)

if nargin < 9, no_flow_segments = []; end
N = size(all_centerlines, 1);
assert(isequal(all_centerlines(:,1), (1:N)'), 'Section IDs must coincide with the row indices 1..N.');
cx = all_centerlines(:,2);
cy = all_centerlines(:,3);
nseg = size(sections_ids, 1);
seg_of_row = zeros(N, 1);
for s = 1:nseg
    seg_of_row(sections_ids(s,1):sections_ids(s,2)) = s;
end

%% Flow fraction and reach type of each segment
bif_child = [];      % first section of every bifurcation branch
bif_sibling = [];    % first section of the other branch of the same bifurcation
for b = 1:size(bif_info, 1)
    ch = reshape(bif_info{b,2}, 1, []);
    bif_child = [bif_child, ch]; %#ok<AGROW>
    bif_sibling = [bif_sibling, fliplr(ch)]; %#ok<AGROW>
end
no_flow_seg = ismember(all_centerlines(sections_ids(:,1), 6), no_flow_segments);
% a branch is a real split only if the other branch of its bifurcation carries discharge
real_split = false(nseg, 1);
for s = 1:nseg
    k = find(bif_child == sections_ids(s,1), 1);
    if ~isempty(k)
        real_split(s) = ~no_flow_seg(seg_of_row(bif_sibling(k)));
    end
end
conf_ids = zeros(1, size(conf_info, 1));
for c = 1:size(conf_info, 1)
    conf_ids(c) = conf_info{c,1};
end
frac = nan(nseg, 1);
frac(network_root) = 1;
for s = reshape(topo_order, 1, [])
    up = predecessor{s};
    if no_flow_seg(s), frac(s) = 0; continue; end
    if s == network_root || isempty(up), continue; end
    fu = frac(up);
    fu = fu(isfinite(fu));
    if isempty(fu), continue; end
    if real_split(s)
        frac(s) = fu(1) / 2;          % branch of a bifurcation
    else
        frac(s) = sum(fu);            % continuation, confluence, or branch whose sibling has no discharge
    end
end
class_seg = repmat({'other multi-thread'}, nseg, 1);
for s = 1:nseg
    if no_flow_seg(s)
        class_seg{s} = 'no flow';
    elseif abs(frac(s) - 1) < 1e-9
        class_seg{s} = 'single channel';
    elseif real_split(s)
        class_seg{s} = 'below bifurcation';
    elseif ismember(sections_ids(s,1), conf_ids)
        class_seg{s} = 'below confluence';
    end
end

%% Stations
st_all = Tpts.Station;
if ~iscell(st_all), st_all = cellstr(st_all); end
X = Tpts.X;
Y = Tpts.Y;
keep = logical(Tpts.valid) & logical(Tpts.inside_corrected_area) & isfinite(X) & isfinite(Y);
stations = unique(st_all(keep));
ns = numel(stations);

nearest_row = zeros(numel(X), 1);
for p = find(keep)'
    [~, nearest_row(p)] = min((cx - X(p)).^2 + (cy - Y(p)).^2);
end

n_points = zeros(ns,1); Xc = nan(ns,1); Yc = nan(ns,1);
seg_main = zeros(ns,1); seg_code = nan(ns,1); row_ref = zeros(ns,1);
share_main = nan(ns,1); dist_axis = nan(ns,1); dist_other = nan(ns,1);
s_from_node = nan(ns,1); flow_fraction = nan(ns,1); reach_class = cell(ns,1);
ord = zeros(ns, 2);
for k = 1:ns
    sel = find(keep & strcmp(st_all, stations{k}));
    segs = seg_of_row(nearest_row(sel));
    useg = unique(segs);
    cnt = arrayfun(@(u) nnz(segs == u), useg);
    [cmax, im] = max(cnt);
    s_main = useg(im);
    xc = mean(X(sel));  yc = mean(Y(sel));
    r_seg = (sections_ids(s_main,1):sections_ids(s_main,2))';
    [d2, kk] = min((cx(r_seg) - xc).^2 + (cy(r_seg) - yc).^2);
    others = seg_of_row ~= s_main;
    r0 = sections_ids(s_main,1);
    rr = r0:r_seg(kk);

    n_points(k)      = numel(sel);
    Xc(k) = xc;  Yc(k) = yc;
    seg_main(k)      = s_main;
    seg_code(k)      = all_centerlines(r_seg(kk), 6);
    row_ref(k)       = r_seg(kk);
    share_main(k)    = cmax / numel(segs);
    dist_axis(k)     = sqrt(d2);
    dist_other(k)    = sqrt(min((cx(others) - xc).^2 + (cy(others) - yc).^2));
    s_from_node(k)   = sum(hypot(diff(cx(rr)), diff(cy(rr))));
    flow_fraction(k) = frac(s_main);
    reach_class{k}   = class_seg{s_main};
    ord(k,:) = [find(topo_order == s_main, 1), r_seg(kk)];
end

[~, io] = sortrows(ord);
ambiguous = share_main < 1;
S = table(stations(io), n_points(io), Xc(io), Yc(io), seg_main(io), seg_code(io), row_ref(io), ...
    share_main(io), ambiguous(io), dist_axis(io), dist_other(io), s_from_node(io), ...
    flow_fraction(io), reach_class(io), ...
    'VariableNames', {'Station','n_points','X','Y','segment','segment_code','section_ID', ...
    'share_main','ambiguous','dist_axis','dist_other','s_from_node','flow_fraction','reach_class'});
end
