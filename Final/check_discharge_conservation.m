function [Tnode, Tseg] = check_discharge_conservation(Q_hm, Q_mm, z_hm, z_mm, all_centerlines, sections_ids, predecessor, bif_info, conf_info, network_root, out_dir)
% CHECK_DISCHARGE_CONSERVATION  Water balance of the flight discharge in the
% hydraulic (HM) and morphodynamic (MM) modules, node by node.
%
% Checks
%   bifurcation   Q(mother, last section) = Q(branch 1, first) + Q(branch 2, first)
%   confluence    Q(first section downstream) = sum of Q(last section) of the upstream segments
%   continuation  Q(first section of a segment) = Q(last section of its only upstream segment)
%   inside a segment, Q must be constant
%   dead ends     segments without a downstream segment: their discharge leaves
%                 the network (correct only for the outlet of the network)
% A node is flagged when the imbalance exceeds 0.5% of the incoming discharge.
%
% INPUTS
%   Q_hm, Q_mm    Q_map of the two modules (one value per section) [m^3/s]
%   z_hm, z_mm    bed elevation of the two modules (z_rect_all, z_rect_morpho)
%   all_centerlines, sections_ids, predecessor, bif_info, conf_info, network_root
%                 network of the hydraulic module (HM_workspace)
%   out_dir       output folder of the csv files
%
% OUTPUTS
%   Tnode  one row per node (bifurcation, confluence, continuation)
%   Tseg   one row per segment: Q at the first and last section of both
%          modules and mean bed difference MM - HM along the segment

Q_hm = Q_hm(:);  Q_mm = Q_mm(:);  z_hm = z_hm(:);  z_mm = z_mm(:);
nseg = size(sections_ids, 1);
code = @(s) all_centerlines(sections_ids(s,1), 6);
seg_of_id = zeros(size(all_centerlines, 1), 1);
for s = 1:nseg
    seg_of_id(sections_ids(s,1):sections_ids(s,2)) = s;
end
tol = 0.005;

%% Segments
first = sections_ids(:,1);  last = sections_ids(:,2);
Qvar_hm = zeros(nseg,1);  Qvar_mm = zeros(nseg,1);  dz = nan(nseg,1);
has_down = false(nseg,1);
for s = 1:nseg
    r = first(s):last(s);
    Qvar_hm(s) = max(Q_hm(r)) - min(Q_hm(r));
    Qvar_mm(s) = max(Q_mm(r)) - min(Q_mm(r));
    dz(s) = mean(z_mm(r) - z_hm(r), 'omitnan');
    for u = reshape(predecessor{s}, 1, [])
        has_down(u) = true;
    end
end
codes = arrayfun(code, (1:nseg)');
Tseg = table(codes, first, last, Q_hm(first), Q_hm(last), Q_mm(first), Q_mm(last), ...
    Qvar_hm, Qvar_mm, 100 * dz, ...
    'VariableNames', {'segment_code','first_ID','last_ID','Q_HM_first','Q_HM_last', ...
    'Q_MM_first','Q_MM_last','Q_HM_range_in_segment','Q_MM_range_in_segment','dz_MM_minus_HM_cm'});

%% Nodes
type = {};  node = [];  up_codes = {};  down_codes = {};
Qin = [];  Qout = [];   % columns: [HM MM]
bif_child = [];
for b = 1:size(bif_info, 1)
    m = bif_info{b,1};  ch = reshape(bif_info{b,2}, 1, []);
    bif_child = [bif_child, ch]; %#ok<AGROW>
    type{end+1,1} = 'bifurcation';  node(end+1,1) = m; %#ok<AGROW>
    up_codes{end+1,1} = sprintf('%d', code(seg_of_id(m))); %#ok<AGROW>
    down_codes{end+1,1} = sprintf('%d ', arrayfun(@(i) code(seg_of_id(i)), ch)); %#ok<AGROW>
    Qin(end+1,:)  = [Q_hm(m), Q_mm(m)]; %#ok<AGROW>
    Qout(end+1,:) = [sum(Q_hm(ch)), sum(Q_mm(ch))]; %#ok<AGROW>
end
for c = 1:size(conf_info, 1)
    d = conf_info{c,1};  up = reshape(conf_info{c,2}, 1, []);
    type{end+1,1} = 'confluence';  node(end+1,1) = d; %#ok<AGROW>
    up_codes{end+1,1} = sprintf('%d ', arrayfun(@(i) code(seg_of_id(i)), up)); %#ok<AGROW>
    down_codes{end+1,1} = sprintf('%d', code(seg_of_id(d))); %#ok<AGROW>
    Qin(end+1,:)  = [sum(Q_hm(up)), sum(Q_mm(up))]; %#ok<AGROW>
    Qout(end+1,:) = [Q_hm(d), Q_mm(d)]; %#ok<AGROW>
end
conf_ids = [];
if ~isempty(conf_info), conf_ids = cell2mat(conf_info(:,1)); end
for s = 1:nseg
    up = predecessor{s};
    if s == network_root || numel(up) ~= 1 || ismember(first(s), bif_child) || ismember(first(s), conf_ids)
        continue
    end
    type{end+1,1} = 'continuation';  node(end+1,1) = first(s); %#ok<AGROW>
    up_codes{end+1,1} = sprintf('%d', code(up)); %#ok<AGROW>
    down_codes{end+1,1} = sprintf('%d', code(s)); %#ok<AGROW>
    Qin(end+1,:)  = [Q_hm(last(up)), Q_mm(last(up))]; %#ok<AGROW>
    Qout(end+1,:) = [Q_hm(first(s)), Q_mm(first(s))]; %#ok<AGROW>
end
imb = Qout - Qin;
flag_hm = abs(imb(:,1)) > tol * abs(Qin(:,1));
flag_mm = abs(imb(:,2)) > tol * abs(Qin(:,2));
Tnode = table(type, node, up_codes, down_codes, Qin(:,1), Qout(:,1), imb(:,1), flag_hm, ...
    Qin(:,2), Qout(:,2), imb(:,2), flag_mm, ...
    'VariableNames', {'type','node_ID','upstream_segments','downstream_segments', ...
    'Q_in_HM','Q_out_HM','imbalance_HM','not_conserved_HM', ...
    'Q_in_MM','Q_out_MM','imbalance_MM','not_conserved_MM'});

%% Report
fprintf('\nWater balance of the flight discharge (imbalance > %.1f%% flagged)\n', 100 * tol);
fprintf('%-13s %6s %10s %10s | %8s %8s %8s | %8s %8s %8s\n', 'node', 'ID', 'up', 'down', ...
    'Qin HM', 'Qout HM', 'diff', 'Qin MM', 'Qout MM', 'diff');
for k = 1:numel(node)
    tag = '';
    if flag_hm(k), tag = [tag '  <-- HM not conserved']; end %#ok<AGROW>
    if flag_mm(k), tag = [tag '  <-- MM not conserved']; end %#ok<AGROW>
    fprintf('%-13s %6d %10s %10s | %8.3f %8.3f %+8.3f | %8.3f %8.3f %+8.3f%s\n', type{k}, node(k), ...
        strtrim(up_codes{k}), strtrim(down_codes{k}), Qin(k,1), Qout(k,1), imb(k,1), ...
        Qin(k,2), Qout(k,2), imb(k,2), tag);
end
bad = find(Qvar_hm > tol * abs(Q_hm(first)) | Qvar_mm > tol * abs(Q_mm(first)));
for s = reshape(bad, 1, [])
    fprintf('Segment %d: Q not constant along the segment (range HM %.3f, MM %.3f m^3/s)\n', ...
        codes(s), Qvar_hm(s), Qvar_mm(s));
end
dead = find(~has_down);
fprintf('\nSegments without a downstream segment (their discharge leaves the network):\n');
for s = reshape(dead, 1, [])
    fprintf('  segment %3d (sections %d-%d): Q HM = %.3f, Q MM = %.3f m^3/s\n', ...
        codes(s), first(s), last(s), Q_hm(last(s)), Q_mm(last(s)));
end
fprintf(['Only one of them should be the outlet of the network: any other one is a branch whose\n' ...
         'reconnection was not found (missing confluence) or a channel that really ends.\n']);

if ~isfolder(out_dir), mkdir(out_dir); end
writetable(Tnode, fullfile(out_dir, 'discharge_balance_nodes.csv'));
writetable(Tseg,  fullfile(out_dir, 'discharge_balance_segments.csv'));
end
