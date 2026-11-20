%% PREPROCESSING & HYDRAULIC MODULE (ORCO)
% Builds the river network from the masks and DTMs of the river segments
% (centerlines, cross-sections, widths, slopes, topology) and computes the
% bed elevation of every cross-section with uniform flow. Segments without
% discharge at the flight (no_flow_segments in the config) get Q = 0 and no
% bed, and their wet pixels keep the raw DTM. Reads every input
% from Input_folder and writes every output to Output_folder:
%   Output_folder/01_Preprocessing     centerlines, widths, topology, figures
%   Output_folder/02_HydraulicModule   riverbed elevation, depth, corrected DTM, DoD
%   Output_folder/04_Validation        GNSS validation of the hydraulic DTM
% All the parameters (and the answers that were typed with input()) are in
% Input_folder/config_ORCO_<year>.m.
%
% Run from MATLAB:   Preprocessing_HydraulicModule
% or from a shell:   matlab -batch "Preprocessing_HydraulicModule"
clear;
clc;
close all;

%% 0. Folders and configuration
root_dir   = fileparts(mfilename('fullpath'));
input_dir  = fullfile(root_dir, 'Input_folder');
output_dir = fullfile(root_dir, 'Output_folder');
config_file = fullfile(input_dir, 'config_ORCO_2025.m');

addpath(root_dir);   % external functions (centerline_from_mask, calculate_widths, ...)
run(config_file);

out_pre   = fullfile(output_dir, '01_Preprocessing');
out_hydro = fullfile(output_dir, '02_HydraulicModule');
out_val   = fullfile(output_dir, '04_Validation', 'HydraulicModule');
out_fig_pre   = fullfile(out_pre, 'figures');
out_fig_hydro = fullfile(out_hydro, 'figures');
for d = {out_pre, out_hydro, out_val, out_fig_pre, out_fig_hydro}
    if ~isfolder(d{1}), mkdir(d{1}); end
end

dtm_files = dir(fullfile(input_dir, dtm_segments_pattern));
mask_files = dir(fullfile(input_dir, mask_segments_pattern));
if isempty(dtm_files) || isempty(mask_files)
    error('No DTM/mask of the river segments found in %s.', input_dir);
end
if numel(dtm_files) ~= numel(mask_files)
    error('%d DTM and %d masks of river segments: they must be paired.', numel(dtm_files), numel(mask_files));
end

%% 1. PRE-PROCESSING: centerlines and cross-sections of each river segment
% Data structures
all_centerlines = []; % cross-sections IDs | 3 centerline coordinates | cross-sections widths | segment code
sections_ids = []; % Start and end IDs of each river segment
sections_ids_with_coords = []; % Start and end IDs with coordinates of each river segment
id_counter = 1;

x_matrix = [];
y_matrix = [];
z_matrix = [];
widths_all = [];

% Iterations for the number of segments
for i = 1:length(dtm_files)

    dtm_path = fullfile(dtm_files(i).folder, dtm_files(i).name);
    mask_path = fullfile(mask_files(i).folder, mask_files(i).name);

    [dtm, R_dtm] = readgeoraster(dtm_path);
    [river_mask, R_mask] = readgeoraster(mask_path);

        dtm = double(dtm);
        % Nodata of the segment DTM (pixels outside the clipping polygon) to NaN:
        % otherwise the bilinear interpolation in calculate_widths mixes the
        % real elevation with -32767 near the edge of the clip and gives false
        % but positive section elevations (e.g. 211.9 m instead of 255.5 m)
        dtm(dtm == -32767) = NaN;

        % Nominal width conversion m --> pixels
    px_size_x = R_mask.CellExtentInWorldX;
    px_size_y = R_mask.CellExtentInWorldY;
    if abs(px_size_x - px_size_y) > 1e-6
        warning('Non-square pixels in mask %s (dx=%.3f, dy=%.3f m): using the mean.', ...
            mask_files(i).name, px_size_x, px_size_y);
    end
    px_size = mean([px_size_x, px_size_y]);   % m/pixel

    Wn_px = Wn / px_size;                     % Nominal width in PIXELS

    if Wn_px < 3
        warning('Wn_px = %.1f px for %s: very low skeleton cleaning threshold, check the resolution.', ...
            Wn_px, mask_files(i).name);
    end

    % River segment edges
    [row, col] = find(river_mask > 0);
    x_coords = R_mask.XWorldLimits(1) + col * R_mask.CellExtentInWorldX;
    y_coords = R_mask.YWorldLimits(2) - row * R_mask.CellExtentInWorldY;

    % River segment orientation
    delta_x = max(x_coords) - min(x_coords);
    delta_y = max(y_coords) - min(y_coords);

    if delta_y > delta_x
        es = 'NS';
    else
        es = 'EW';
    end

    % Centerline calculation
    [cl, ~] = centerline_from_mask(river_mask, es, Wn_px, plotornot);

    [x_cl, y_cl] = intrinsicToWorld(R_mask, cl(:, 1), cl(:, 2));
    x_world = linspace(R_mask.XWorldLimits(1), R_mask.XWorldLimits(2), size(river_mask, 2));
    y_world = linspace(R_mask.YWorldLimits(2), R_mask.YWorldLimits(1), size(river_mask, 1));


    z_cl = interp2(x_world, y_world, dtm, x_cl, y_cl, 'linear');
    xyz_cl = [x_cl, y_cl, z_cl];

    dist_cum = [0; cumsum(sqrt(sum(diff(xyz_cl(:, 1:2)).^2, 2)))];
    dist_interp = 0:dx:dist_cum(end);
    x_interp = interp1(dist_cum, xyz_cl(:, 1), dist_interp, 'linear');
    y_interp = interp1(dist_cum, xyz_cl(:, 2), dist_interp, 'linear');
    z_interp = interp1(dist_cum, xyz_cl(:, 3), dist_interp, 'linear');
    for j = 1:length(z_interp)
        if z_interp(j) < 0
            if j > 1 && z_interp(j-1) >= 0
                z_interp(j) = z_interp(j-1);
            elseif j < length(z_interp) && z_interp(j+1) >= 0
                z_interp(j) = z_interp(j+1);  % otherwise take the next value
            end
        end
    end

    % Verify upstream-downstream direction of river segments
    if z_interp(1) < z_interp(end)
        x_interp = flipud(x_interp);
        y_interp = flipud(y_interp);
        z_interp = flipud(z_interp);
    end

    num_points = length(x_interp);
    [~, mask_name, ~] = fileparts(mask_files(i).name);
    mask_code = extractAfter(mask_name, 'I_'); % Cross-section number extraction (es. '000','001', etc)
    segment_code = str2double(mask_code); % (000 -> 0, etc.)
    segment_codes = repmat(segment_code, num_points, 1); % same code for all segment's points
    % Assign IDs
    ids = (id_counter:id_counter + num_points - 1)';
    id_counter = id_counter + num_points;

    % Width calculation
    [widths, sect_x, sect_y, sect_z, z_sect_mean] = calculate_widths(...
        x_interp, y_interp, z_interp, x_world, y_world, river_mask, dtm);
    widths_all = vertcat(widths_all, widths(:));
    x_matrix = vertcat(x_matrix, sect_x);
    y_matrix = vertcat(y_matrix, sect_y);
    z_matrix = vertcat(z_matrix, sect_z);

    % Plot binary mask of each river segment with centerline and cross-sections
    fig = figure;
    imagesc(R_mask.XWorldLimits, R_mask.YWorldLimits, flipud(river_mask));
    colormap('gray');
    set(gca, 'YDir', 'normal');
    hold on;

    % Plot cross-sections
    [r, ~] = size(sect_x);
    for j = 1:r
    plot(sect_x(j, :), sect_y(j, :), '.r');
    hold on;
    end

    % Plot centerline
    plot(x_interp, y_interp, 'ob');

    title(sprintf('River segment %s: Centerline and cross-sections', mask_code));
    xlabel('X');
    ylabel('Y');
    axis equal;
    hold off;
    drawnow;
    save_fig(fig, save_figures, out_fig_pre, sprintf('PRE_segment_%s_centerline_sections', mask_code));

    % Be sure to have column vectors
    ids = ids(:);
    x_interp = x_interp(:);
    y_interp = y_interp(:);
    z_interp = z_interp(:);
    widths = widths(:);

    widths(isnan(widths)) = 0;
    all_centerlines = [all_centerlines; ...
    ids, x_interp, y_interp, z_sect_mean, widths, segment_codes];
    sections_ids = [sections_ids; ids(1), ids(end)];
    start_coords = [x_interp(1), y_interp(1), z_interp(1)];
    end_coords = [x_interp(end), y_interp(end), z_interp(end)];
    sections_ids_with_coords = [sections_ids_with_coords; ...
        ids(1), start_coords, ids(end), end_coords];
end

disp('Done!');

%% 2D Plot
x_coords = all_centerlines(:, 2);
y_coords = all_centerlines(:, 3);
z_coords = all_centerlines(:, 4);

x_start = sections_ids_with_coords(:, 2);
y_start = sections_ids_with_coords(:, 3);
x_end = sections_ids_with_coords(:, 6);
y_end = sections_ids_with_coords(:, 7);

fig = figure;
scatter(x_coords, y_coords, 10, 'b', 'filled');
title('Centerlines points');
xlabel('X');
ylabel('Y');
axis equal;
grid on;
hold on;
scatter(x_start, y_start, 10, 'or');
scatter(x_end, y_end, 10, 'og');
legend('centerline', 'start', 'end');

    mask_codes = strings(length(mask_files), 1);

for i = 1:length(mask_files)
    [~, mask_name, ~] = fileparts(mask_files(i).name);
    mask_codes(i) = extractAfter(mask_name, 'I_');
end

if length(mask_codes) ~= length(x_start)
    error('Codes number is not equal to river segments!');
end

for i = 1:length(x_start)
    text(x_start(i), y_start(i), mask_codes(i), ...
    'FontSize', 8, 'Color', 'r', ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right');
end
hold off;
save_fig(fig, save_figures, out_fig_pre, 'PRE_centerlines_before_flip');

%% Flip river segments UPSTREAM - DOWNSTREAM
% Codes from config (segments_to_flip); interactive = true asks them with input()
if interactive
    invert_input = input('Are there river segments to flip (yes/no): ', 's');
    if strcmpi(invert_input, 'yes')
        tratti_da_invertire = input('Write codes of river segments to flip (spaced): ', 's');
        tratti_da_invertire = strsplit(tratti_da_invertire);
    else
        tratti_da_invertire = {};
    end
else
    tratti_da_invertire = segments_to_flip;
end

for i = 1:length(tratti_da_invertire)
    codice_tratto = tratti_da_invertire{i};

    % Find river segment ID
    idx = find(mask_codes == codice_tratto);

    if ~isempty(idx)
        % Flip
        start_idx = find(all_centerlines(:,1) == sections_ids(idx,1));
        end_idx = find(all_centerlines(:,1) == sections_ids(idx,2));
        range = min(start_idx, end_idx):max(start_idx, end_idx);
        all_centerlines(range, 2:end) = flipud(all_centerlines(range, 2:end));

        % Flip ALSO the cross-section matrices and the widths, otherwise they
        % keep the pre-flip order and are misaligned with all_centerlines
        x_matrix(range, :)  = flipud(x_matrix(range, :));
        y_matrix(range, :)  = flipud(y_matrix(range, :));
        z_matrix(range, :)  = flipud(z_matrix(range, :));
        widths_all(range)   = flipud(widths_all(range));

        % Update coordinates
        sections_ids_with_coords(idx, [2:4, 6:8]) = sections_ids_with_coords(idx, [6:8, 2:4]);
    else
        warning('River segment %s to flip not found among the masks.', codice_tratto);
    end
end
x_coords = all_centerlines(:, 2);
y_coords = all_centerlines(:, 3);
z_coords = all_centerlines(:, 4);
x_start = sections_ids_with_coords(:, 2);
y_start = sections_ids_with_coords(:, 3);
x_end = sections_ids_with_coords(:, 6);
y_end = sections_ids_with_coords(:, 7);
self_dist = nan(size(all_centerlines,1),1);
for row = 1:size(all_centerlines,1)
    dx = x_matrix(row,:) - all_centerlines(row,2);
    dy = y_matrix(row,:) - all_centerlines(row,3);
    self_dist(row) = min(sqrt(dx.^2 + dy.^2), [], 'omitnan');
end
n_bad = nnz(self_dist > 5);
if n_bad > 0
    error('Centerline/sections misalignment for %d rows after flip: verify!', n_bad);
end

%% Same plot to verify
fig = figure;
scatter(x_coords, y_coords, 10, 'b', 'filled');
title('Centerlines points');
xlabel('X');
ylabel('Y');
axis equal;
grid on;
hold on;
scatter(x_start, y_start, 10, 'or');
scatter(x_end, y_end, 10, 'og');
legend('centerline', 'start', 'end');

for i = 1:length(x_start)
    text(x_start(i), y_start(i), mask_codes(i), ...
    'FontSize', 8, 'Color', 'r', ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right');
end
hold off;
save_fig(fig, save_figures, out_fig_pre, 'PRE_centerlines_after_flip');

%% ADJACENCY MATRIX
% River segments number
num_sections = size(sections_ids, 1);

adjacency_matrix = zeros(num_sections * 2 + 1);

all_ids = [sections_ids(:, 1); sections_ids(:, 2)]; % start and end IDs
adjacency_matrix(2:end, 1) = all_ids;
adjacency_matrix(1, 2:end) = all_ids;

% Loop to fill adjacency matrix (distance_threshold from config)
for i = 1:num_sections
    start_coords = [x_start(i), y_start(i)]; % [x, y] start river segment i

    for j = 1:num_sections
            if i == j
            continue;
        end

        end_coords = [x_end(j), y_end(j)]; % [x, y] end of river segment j
        dist = sqrt(sum((start_coords - end_coords).^2)); % Euclidean distance

        if dist < distance_threshold
            % Find corrisponding IDs in the adjacency matrix
            row_i = find(adjacency_matrix(2:end, 1) == sections_ids(j, 2)) + 1; % Row of the ending point j
            col_j = find(adjacency_matrix(1, 2:end) == sections_ids(i, 1)) + 1; % Column of the starting point i

            % Update adjacency matrix
            adjacency_matrix(row_i, col_j) = 1;  % The ending point j is close to the starting point i
            adjacency_matrix(col_j, row_i) = -1; % The starting point i is close to the ending point j
        end
    end
end

%% Bifurcations and confluences from the adjacency matrix
% Cell array initialization to store bif and conf details
bif_info = {};  % Each row: {SectionID, [Downstream IDs]}
conf_info = {}; % Each row: {SectionID, [Upstream IDs]}
bif_info_tratti = {};   % {SegmentCode, [Downstream SegmentCodes]}
conf_info_tratti = {};  % {SegmentCode, [Upstream SegmentCodes]}

num_ids = size(adjacency_matrix, 1) - 1; % exclude row and column title

% Column analysis (from 2nd to last)
for col = 2:(num_ids+1)

    sectionID = adjacency_matrix(1, col); % section ID

    % Find 1
    % SectionID gets from more than one section --> confluence
    contributorIndices = find(adjacency_matrix(2:end, col) == 1);
    if length(contributorIndices) > 1
        contributorIDs = adjacency_matrix(contributorIndices + 1, 1);
        conf_info{end+1, 1} = sectionID;
        conf_info{end, 2} = contributorIDs;
        % Find river segments codes
        code_current = all_centerlines(all_centerlines(:,1) == sectionID, 6);
        code_contributors = zeros(length(contributorIDs), 1);
        for k = 1:length(contributorIDs)
            idx = find(all_centerlines(:,1) == contributorIDs(k), 1);
            code_contributors(k) = all_centerlines(idx, 6);
        end
        conf_info_tratti{end+1, 1} = code_current;
        conf_info_tratti{end, 2} = code_contributors;
    end

    % Find -1
    % SectionID feeds more than one section --> bifurcation
    feederIndices = find(adjacency_matrix(2:end, col) == -1); % rows that this section flows into
    if length(feederIndices) > 1
        feederIDs = adjacency_matrix(feederIndices + 1, 1);
        bif_info{end+1, 1} = sectionID;
        bif_info{end, 2} = feederIDs;
        % Find river segments codes
        code_current = all_centerlines(all_centerlines(:,1) == sectionID, 6);
        code_feeders = zeros(length(feederIDs), 1);
        for k = 1:length(feederIDs)
            idx = find(all_centerlines(:,1) == feederIDs(k), 1);
            code_feeders(k) = all_centerlines(idx, 6);
        end
        bif_info_tratti{end+1, 1} = code_current;
        bif_info_tratti{end, 2} = code_feeders;
    end
end

%% Overcome bridges issue-manually insert bifurcation and confluences not found
% From config (forced_bifurcations, forced_confluences); interactive = true
% asks them with input()
if interactive
    forced_bifurcations = zeros(0, 3);
    while true
        str = lower(strtrim(input('Are there other bifurcations? (yes/no): ', 's')));
        if strcmp(str, 'no')
            break;
        elseif strcmp(str, 'yes')
            tokens = strsplit(input('Write river segments codes spaced (bifurcating_branch branch1 branch2): ', 's'));
            if length(tokens) ~= 3
                warning('Invalid input. You must specify exactly 3 codes.');
                continue;
            end
            forced_bifurcations(end+1, :) = str2double(tokens); %#ok<SAGROW>
        else
            disp('Please answer with "yes" or "no".');
        end
    end
    forced_confluences = zeros(0, 3);
    while true
        str = lower(strtrim(input('Are there other confluences? (yes/no): ', 's')));
        if strcmp(str, 'no')
            break;
        elseif strcmp(str, 'yes')
            tokens = strsplit(input('Write river segments codes spaced (branch1 branch2 merging_branch): ', 's'));
            if length(tokens) ~= 3
                warning('Invalid input. You must specify exactly 3 codes.');
                continue;
            end
            forced_confluences(end+1, :) = str2double(tokens); %#ok<SAGROW>
        else
            disp('Please answer with "yes" or "no".');
        end
    end
end

for b = 1:size(forced_bifurcations, 1)
    % Parse segment codes
    main_code = forced_bifurcations(b, 1);
    branch1_code = forced_bifurcations(b, 2);
    branch2_code = forced_bifurcations(b, 3);

    % Get end ID of the bifurcation and start IDs of the two branches
    main_idx = find(all_centerlines(:, 6) == main_code);
    b1_idx = find(all_centerlines(:, 6) == branch1_code);
    b2_idx = find(all_centerlines(:, 6) == branch2_code);

    if isempty(main_idx) || isempty(b1_idx) || isempty(b2_idx)
        warning('Could not find all river segments in all_centerlines (bifurcation %s).', mat2str(forced_bifurcations(b, :)));
        continue;
    end

    % Add river segments codes to bif_info_tratti
    bif_info_tratti(end+1, :) = {main_code, [branch1_code, branch2_code]};

    main_last_id = all_centerlines(main_idx(end), 1);
    b1_first_id = all_centerlines(b1_idx(1), 1);
    b2_first_id = all_centerlines(b2_idx(1), 1);

    % Add points IDs to bif_info
    bif_info(end+1, :) = {main_last_id, [b1_first_id, b2_first_id]};
end

for c = 1:size(forced_confluences, 1)
    % Parse segment codes
    branch1_code = forced_confluences(c, 1);
    branch2_code = forced_confluences(c, 2);
    merging_code = forced_confluences(c, 3);

    % Get end IDs of the two branches and start ID of the confluence
    b1_idx = find(all_centerlines(:, 6) == branch1_code);
    b2_idx = find(all_centerlines(:, 6) == branch2_code);
    merge_idx = find(all_centerlines(:, 6) == merging_code);

    if isempty(b1_idx) || isempty(b2_idx) || isempty(merge_idx)
        warning('Could not find all river segments in all_centerlines (confluence %s).', mat2str(forced_confluences(c, :)));
        continue;
    end

    % Add river segments codes to conf_info_tratti
    conf_info_tratti(end+1, :) = {merging_code, [branch1_code, branch2_code]};

    b1_last_id = all_centerlines(b1_idx(end), 1);
    b2_last_id = all_centerlines(b2_idx(end), 1);
    merge_first_id = all_centerlines(merge_idx(1), 1);

    % Add points IDs to conf_info
    conf_info(end+1, :) = {merge_first_id, [b1_last_id, b2_last_id]};
end

%% ============================================================
%% TOPOLOGICAL ORDERING (fix for discharge propagation)
%% ============================================================
n_points = size(all_centerlines, 1);
section_to_index = containers.Map(all_centerlines(:,1), 1:n_points);
num_sections = size(sections_ids, 1);
start_ids = sections_ids(:,1);
end_ids   = sections_ids(:,2);

predecessor = cell(num_sections, 1);   % upstream segment(s) of each segment (2 at a confluence)
edges = [];                            % graph edges [upstream segment, downstream segment]

for s = 1:num_sections
    start_id = start_ids(s);

    % 1) Does the segment start at a CONFLUENCE?
    conf_match = [];
    if ~isempty(conf_info)
        conf_match = find(cellfun(@(x) isequal(x, start_id), conf_info(:,1)));
    end
    if ~isempty(conf_match)
        up_ids = conf_info{conf_match, 2};
        for k = 1:numel(up_ids)
            s_up = find(end_ids == up_ids(k));
            if isempty(s_up)
                error('Topology: upstream segment not found for confluence (ID %d).', up_ids(k));
            end
            predecessor{s} = [predecessor{s}, s_up];
            edges = [edges; s_up, s];
        end
        continue
    end

    % 2) Is the segment a BIFURCATION branch?
    is_branch = false;
    for b = 1:size(bif_info,1)
        if any(bif_info{b,2} == start_id)
            s_up = find(end_ids == bif_info{b,1});
            if isempty(s_up)
                error('Topology: upstream segment not found for bifurcation (ID %d).', bif_info{b,1});
            end
            predecessor{s} = s_up;
            edges = [edges; s_up, s];
            is_branch = true;
            break
        end
    end
    if is_branch, continue, end

    % 3) Simple continuation: single connection in adjacency_matrix
    col = find(adjacency_matrix(1,2:end) == start_id) + 1;
    if ~isempty(col)
        row = find(adjacency_matrix(2:end, col) == 1);
        if numel(row) == 1
            up_id = adjacency_matrix(row+1, 1);
            s_up = find(end_ids == up_id);
            if isempty(s_up)
                error('Topology: upstream segment not found for simple continuation (ID %d).', up_id);
            end
            predecessor{s} = s_up;
            edges = [edges; s_up, s];
            continue
        elseif numel(row) > 1
            error('Topology: segment %d has several upstream connections not registered as a confluence.', s);
        end
    end

    % 4) No connection found -> network head
    predecessor{s} = [];
end

% Topological ordering
if isempty(edges)
    topo_order = 1:num_sections;
else
    G = digraph(edges(:,1), edges(:,2), [], num_sections);
    try
        topo_order = toposort(G);
    catch ME
        error('The segment graph contains a cycle: check bif_info/conf_info. %s', ME.message);
    end
end

% --- Diagnostic check: comparison with the assumption "increasing code = downstream"
codes_by_section = arrayfun(@(s) all_centerlines(find(all_centerlines(:,1)==sections_ids(s,1),1), 6), 1:num_sections);
[~, order_by_code] = sort(codes_by_section);
if ~isequal(order_by_code(:), topo_order(:))
    % Comparison by RANK (position), not by set: setdiff is useless here
    % because the two orderings always contain the same segments.
    rank_by_code = zeros(num_sections,1);
    rank_by_code(order_by_code) = 1:num_sections;
    rank_topo = zeros(num_sections,1);
    rank_topo(topo_order) = 1:num_sections;

    mismatched = find(rank_by_code ~= rank_topo);
    mismatch_codes = codes_by_section(mismatched);

    warning(['The order by increasing code does NOT match the actual topological order. ' ...
             'Using the geometric (correct) order. Segment codes in a different position: %s'], ...
             mat2str(mismatch_codes));
end

% --- Network head(s) and upstream discharge
roots = find(cellfun(@isempty, predecessor));
if isempty(roots)
    error('Topology: no network head found (possible cycle or data error).');
end
if numel(roots) > 1
    warning(['Found %d network heads with no upstream segment: %s. ' ...
             'Check which one is the upstream gauging station (San Benigno) and whether the others ' ...
             'are tributaries that need a known discharge assigned manually.'], ...
             numel(roots), mat2str(roots));
    root_codes = codes_by_section(roots);
    disp(table(roots(:), root_codes(:), 'VariableNames', {'row_index','segment_code'}));
end
network_root = roots(1);

%% Clean isolated glitches of z_DTM (same rule for every segment)
% Replaces the old curvature filter (threshold 1 m) and the manual list
% segments_to_smooth (moving median on selected segments), removed on
% 20/11/2026: with the nodata of the segment DTMs converted to NaN no segment
% oscillates any more, and the moving median only removed a few isolated
% spikes in the listed segments while flattening real downstream steps.
%
% Rule. For every interior section i of a segment:
%   e_i = z_i - (z_{i-1} + z_{i+1})/2       (excess over the neighbours)
% section i is an isolated glitch if
%   1) it is a local extremum: (z_i - z_{i-1})*(z_i - z_{i+1}) >= 0, so a
%      step (z_i between its neighbours) is never touched, and
%   2) |e_i| > max(glitch_tol, 3*sigma_e), with sigma_e = 1.4826*MAD(e) of the
%      segment (robust spread of the excess: in a noisier segment only the
%      outliers are corrected).
% A glitch is replaced by the mean of its two raw neighbours. The test and
% the replacement use the raw values (z_original), so that a correction does
% not propagate to the neighbouring sections.
% With the 19/11/2026 data and glitch_tol = 0.20 m the rule corrects 4
% sections: 14 (segment 0, +0.55 m), 492 (segment 19, -0.38 m), 728
% (segment 26, +0.38 m), 1239 (segment 35, +0.53 m).
z_original = all_centerlines(:,4);   % raw z_DTM, snapshot before the correction
outlier_flags = false(height(all_centerlines), 1);
glitch_excess = nan(height(all_centerlines), 1);
tronchi_uniq = unique(all_centerlines(:,6));
for t = tronchi_uniq'
    idx = find(all_centerlines(:,6) == t);
    z = z_original(idx);
    if numel(z) < 5, continue; end
    e = nan(size(z));
    e(2:end-1) = z(2:end-1) - (z(1:end-2) + z(3:end)) / 2;
    is_ext = false(size(z));
    is_ext(2:end-1) = (z(2:end-1) - z(1:end-2)) .* (z(2:end-1) - z(3:end)) >= 0;
    ei = e(2:end-1);
    sigma_e = 1.4826 * median(abs(ei - median(ei, 'omitnan')), 'omitnan');
    soglia = max(glitch_tol, 3 * sigma_e);
    out = is_ext & abs(e) > soglia;
    outlier_flags(idx(out)) = true;
    glitch_excess(idx(out)) = e(out);
end
fprintf('%d sections with an isolated glitch of z_coords (glitch_tol = %.2f m)\n', nnz(outlier_flags), glitch_tol);
disp(table(all_centerlines(outlier_flags,1), all_centerlines(outlier_flags,6), z_original(outlier_flags), ...
    glitch_excess(outlier_flags), 'VariableNames', {'ID','segment','z_raw','excess'}))

for row = find(outlier_flags)'
    idx = find(all_centerlines(:,6) == all_centerlines(row,6));
    pos_in_tronco = find(idx == row);
    z_prev = z_original(idx(pos_in_tronco-1));
    z_next = z_original(idx(pos_in_tronco+1));
    all_centerlines(row, 4) = mean([z_prev, z_next]);
end
% Re-synchronise z_coords with all_centerlines(:,4) after the glitch cleaning,
% otherwise slopes and hydraulics stay tied to the uncleaned version
z_coords = all_centerlines(:, 4);

%% Slope calculation (mean of segment-by-segment slopes)
slope_sect = zeros(size(all_centerlines, 1), 1);
n_sections = size(sections_ids_with_coords, 1);

for s = 1:n_sections
    idx_start = sections_ids_with_coords(s, 1);
    idx_end   = sections_ids_with_coords(s, 5);

    if (idx_end - idx_start + 1) >= 2 % If the river segment has more than 2 points
        % Extract centerline points in that river segment
        x = all_centerlines(idx_start:idx_end, 2);
        y = all_centerlines(idx_start:idx_end, 3);
        z = all_centerlines(idx_start:idx_end, 4);

        % Calculate differences among consecutive points
        dz = diff(z);
        dx = diff(x);
        dy = diff(y);
        ds = sqrt(dx.^2 + dy.^2); % Euclidean distance
        ss=zeros(1,numel(ds)+1);
        for j=2:numel(ds)+1
            ss(j)=ss(j-1)+ds(j-1);
        end

        p=polyfit(ss,z.',1);
        slope_sect(idx_start:idx_end) =max(slope_min,atan(-p(1)));


    else
        warning("River segment %d is too short. Defined slope assigned", s);
        slope_sect(idx_start:idx_end) = slope_min;
    end
end

sect_slope = slope_sect(:);

%% Save PRE-PROCESSING outputs
% all_centerlines columns: [ID, X, Y, z_DTM, W mask, segment code]
T_cl = array2table([all_centerlines, slope_sect], 'VariableNames', ...
    {'ID','X','Y','z_DTM','W_mask','segment_code','slope'});
writetable(T_cl, fullfile(out_pre, sprintf('PRE_centerlines_%s.csv', run_tag)));
T_seg = array2table([sections_ids_with_coords, codes_by_section(:)], 'VariableNames', ...
    {'ID_start','X_start','Y_start','Z_start','ID_end','X_end','Y_end','Z_end','segment_code'});
writetable(T_seg, fullfile(out_pre, sprintf('PRE_segments_endpoints_%s.csv', run_tag)));
save(fullfile(out_pre, sprintf('PRE_workspace_%s.mat', run_tag)), ...
    'all_centerlines', 'sections_ids', 'sections_ids_with_coords', 'mask_codes', ...
    'x_matrix', 'y_matrix', 'z_matrix', 'widths_all', 'slope_sect', ...
    'adjacency_matrix', 'bif_info', 'conf_info', 'bif_info_tratti', 'conf_info_tratti', ...
    'predecessor', 'topo_order', 'network_root', 'codes_by_section', 'outlier_flags');
disp('Pre-processing outputs saved.');

%% 2. HYDRAULIC MODULE
% g, n_manning, Q_upstream, level_upstream from config
h_rect_all = NaN(n_points,1); % Depth
z_rect_all = NaN(n_points,1); % riverbed elevation matrix
sectionWet(1:n_points) = struct('h', NaN, 'quot', NaN, 'Q', NaN); % structure
Q_map = NaN(n_points, 1);       % Discharge preallocation for each section
root_start_id = sections_ids(network_root, 1);
Q_map(section_to_index(root_start_id)) = Q_upstream;
% IDs of bifurcation/confluence nodes (empty if there are none)
bif_ids  = zeros(0,1);
conf_ids = zeros(0,1);
if ~isempty(bif_info),  bif_ids  = cell2mat(bif_info(:,1));  end
if ~isempty(conf_info), conf_ids = cell2mat(conf_info(:,1)); end
% Sections of the segments without discharge at the flight (no_flow_segments)
no_flow_row = ismember(all_centerlines(:,6), no_flow_segments);
%% For cycle for all sections
row_idx_cell = arrayfun(@(s) find(all_centerlines(:,1) >= sections_ids(s,1) & ...
                     all_centerlines(:,1) <= sections_ids(s,2)), topo_order, 'UniformOutput', false);
iter_rows = vertcat(row_idx_cell{:});
is_segment_start = ismember(all_centerlines(:,1), sections_ids(:,1));

for k = 1:n_points
    i = iter_rows(k);
    % ----------------------- SEGMENT WITHOUT DISCHARGE -----------------------
    % Q = 0 at the flight: depth and bed are not computed (NaN), so the
    % section is left out of the bed interpolation and of the DTM correction.
    if no_flow_row(i)
        Q_map(i) = 0;
        continue
    end
    % ----------------------------- BIFURCATION -----------------------------
    if ismember(all_centerlines(i,1), bif_ids) % If bifurcation
        width = widths_all(i);
        slope = sect_slope(i);
        Q_map(i) = Q_map(i-1);
        h_rect = (Q_map(i) / (width * (1/n_manning) * sqrt(slope)))^(3/5);
        z_rect = z_coords(i) - h_rect;
        sectionWet(i).h = h_rect;
        sectionWet(i).quot = z_rect;
        sectionWet(i).Q = Q_map(i);
        h_rect_all(i) = h_rect;
        z_rect_all(i) = z_rect;

        % Find the two branches of the bifurcation
        bif_idx = find(bif_ids == all_centerlines(i,1));
        next_sections = bif_info{bif_idx, 2};

            if length(next_sections) == 2
                sec1 = next_sections(1);
                sec2 = next_sections(2);

                % Bifurcation branches indexes
                idx1 = section_to_index(sec1);
                idx2 = section_to_index(sec2);
                % (if sect1 is the 5th ID in all_centerlines(:,1) follows idx1=5)

                width1 = widths_all(idx1);
                width2 = widths_all(idx2);
                slope1 = sect_slope(idx1);
                slope2 = sect_slope(idx2);

                % System of eq: [Q1, Q2, h1, h2]
                F = @(x) [
                    % 1. Discharge conservation
                    x(1) + x(2) - Q_map(i);
                    % 2. Chezy's eq 1
                    x(1) - (1/n_manning) * width1 * x(3)^(5/3) * sqrt(slope1);
                    % 3. Chezy's eq 2
                    x(2) - (1/n_manning) * width2 * x(4)^(5/3) * sqrt(slope2);
                    % 4. Bernoulli
                    x(3) + x(1)^2 / (2 * g * (width1 * x(3))^2) ...
                          - (x(4) + x(2)^2 / (2 * g * (width2 * x(4))^2));
                ];

                no_flow_child = no_flow_row([idx1 idx2]);
                if any(no_flow_child)
                    % One branch carries no discharge at the flight: the whole
                    % discharge goes to the other branch (uniform flow), no
                    % node system. The branch without discharge gets Q = 0 and
                    % no depth (NaN).
                    Q1 = Q_map(i) * double(~no_flow_child(1));
                    Q2 = Q_map(i) - Q1;
                    h_rect1 = NaN;  h_rect2 = NaN;
                    if Q1 > 0, h_rect1 = (Q1 / (width1 * (1/n_manning) * sqrt(slope1)))^(3/5); end
                    if Q2 > 0, h_rect2 = (Q2 / (width2 * (1/n_manning) * sqrt(slope2)))^(3/5); end
                    fprintf('Bifurcation %d: branch without discharge (no_flow_segments), whole discharge to the other branch\n', ...
                        all_centerlines(i,1));
                else
                    idx_current = section_to_index(all_centerlines(i,1));
                    if ~all(isnan(h_rect_all(idx_current)))
                        h_init = [h_rect_all(idx_current); h_rect_all(idx_current)];
                    else
                        h_init = [level_upstream; level_upstream];
                    end

                    x0 = [Q_map(i)/2; Q_map(i)/2; h_init(1); h_init(2)];

                    options = optimoptions('fsolve', 'Display', 'iter', 'FunctionTolerance', 1e-8);

                    [x_sol, fval, exitflag] = fsolve(F, x0, options);

                    if exitflag > 0 && all(isreal(x_sol)) && all(x_sol > 0)
                        Q1 = x_sol(1); Q2 = x_sol(2);
                        h_rect1 = x_sol(3); h_rect2 = x_sol(4);
                    else
                        warning('No physical solution for bifurcation in section %d', all_centerlines(i,1));
                        disp('------ BIFURCATION FAILURE INFO ------');
                        disp(['Section ID: ', num2str(all_centerlines(i,1))]);
                        disp(['Widths: ', num2str([width1, width2])]);
                        disp(['Slopes: ', num2str([slope1, slope2])]);
                        disp(['Initial guess: ', num2str(x0')]);
                        disp(['Solution: ', num2str(x_sol')]);
                        disp(['fval: ', num2str(fval')]);
                        % Fallback:one of the two branches is almost flat
                        Q1 = Q_map(i) * width1 / (width1 + width2);
                        Q2 = Q_map(i) - Q1;

                        h_rect1 = (Q1 / (width1 * (1/n_manning) * sqrt(slope1)))^(3/5);
                        h_rect2 = (Q2 / (width2 * (1/n_manning) * sqrt(slope2)))^(3/5);
                    end
                end   % branch without discharge / node system
                % Save
                    Q_map(idx1) = Q1;
                    Q_map(idx2) = Q2;
                    z_rect1 = z_coords(idx1) - h_rect1;
                    z_rect2 = z_coords(idx2) - h_rect2;

                    sectionWet(idx1).h = h_rect1;
                    sectionWet(idx1).quot = z_rect1;
                    sectionWet(idx1).Q = Q1;
                    sectionWet(idx2).h = h_rect2;
                    sectionWet(idx2).quot = z_rect2;
                    sectionWet(idx2).Q = Q2;
                    h_rect_all(idx1) = h_rect1;
                    h_rect_all(idx2) = h_rect2;
                    z_rect_all(idx1) = z_rect1;
                    z_rect_all(idx2) = z_rect2;
            end
    % ----------------------------- CONFLUENCE -----------------------------
    elseif ismember(all_centerlines(i,1), conf_ids) % If confluence
        width = widths_all(i);
        slope = sect_slope(i);

        conf_idx = find(conf_ids == all_centerlines(i,1));
        upstream_sections = conf_info{conf_idx, 2};

        idx1 = section_to_index(upstream_sections(1));
        idx2 = section_to_index(upstream_sections(2));

        Q1 = Q_map(idx1);
        Q2 = Q_map(idx2);
        Q_total = Q1 + Q2;
        Q_map(i) = Q_total;
        h_rect = (Q_total / (width * (1/n_manning) * sqrt(slope)))^(3/5);
        z_rect = z_coords(i) - h_rect;
        sectionWet(i).h = h_rect;
        sectionWet(i).quot = z_rect;
        sectionWet(i).Q = Q_total;
        h_rect_all(i) = h_rect;
        z_rect_all(i) = z_rect;
        % propagate only inside the same segment: the first section of another
        % segment takes its discharge from its own upstream segment
        if i < n_points && ...
           all_centerlines(i+1,6) == all_centerlines(i,6) && ...
           ~ismember(all_centerlines(i+1,1), conf_ids) && ...
           ~ismember(all_centerlines(i+1,1), bif_ids)
            Q_map(i+1) = Q_total;
        end
    else
% ----------------------------- SINGLE CHANNEL -----------------------------
        width = widths_all(i);
        slope = sect_slope(i);
        if isnan(width) || width <= 0 || isnan(slope) || slope <= 0
           warning('Invalid width or slope at section %d', i);
           continue;
        end
        if isnan(Q_map(i))
            if is_segment_start(i)
                seg_idx = find(sections_ids(:,1) == all_centerlines(i,1));
                up = predecessor{seg_idx};
                if isempty(up)
                    warning('Segment %d (network head): upstream discharge not defined. Section skipped.', seg_idx);
                    continue
                end
                up_row = section_to_index(sections_ids(up,2));
                if isnan(Q_map(up_row))
                    error('Segment %d: discharge of the upstream segment (%d) not computed yet. Topological order error.', seg_idx, up);
                end
                Q_map(i) = Q_map(up_row);
            elseif k > 1
                Q_map(i) = Q_map(iter_rows(k-1));   % continuation inside the same segment
            else
                warning('Q_map(%d) not defined and the previous value cannot be retrieved.', i);
                continue;
            end
        end

        Q_local = Q_map(i);
        h_rect = (Q_local / (width * (1/n_manning) * sqrt(slope)))^(3/5);
        z_rect = z_coords(i) - h_rect;
        sectionWet(i).h = h_rect;
        sectionWet(i).quot = z_rect;
        sectionWet(i).Q = Q_local;
        h_rect_all(i) = h_rect;
        z_rect_all(i) = z_rect;
        % Discharge propagation, only inside the same segment: the first
        % section of another segment takes its discharge from its own upstream
        % segment (otherwise it would inherit the discharge of the previous
        % row, i.e. of an unrelated segment)
        if i < n_points && ...
           all_centerlines(i+1,6) == all_centerlines(i,6) && ...
           ~ismember(all_centerlines(i+1,1), conf_ids) && ...
           ~ismember(all_centerlines(i+1,1), bif_ids)

           Q_map(i+1) = Q_local;
        end
    end
end

%% Save HYDRAULIC MODULE outputs (centerline)
% Save centerline with corrected bathymetry in CSV format
fileID = fopen(fullfile(out_hydro, sprintf('HM_centerline_riverbed_%s.csv', run_tag)), 'w');

fprintf(fileID, 'X,Y,Z\n');

for i = 1:length(x_coords)
    fprintf(fileID, '%.3f,%.3f,%.3f\n', x_coords(i), y_coords(i), z_rect_all(i));
end

fclose(fileID);
% Save depth for each cross-section
save(fullfile(out_hydro, sprintf('HM_h_rect_all_%s.mat', run_tag)), 'h_rect_all');
% Create a bathymetry corrected matrix like z_matrix
z_rect_all_matrix = NaN(size(z_matrix));

for i = 1:size(z_matrix, 1)
    mask = ~isnan(z_matrix(i, :));
    z_rect_all_matrix(i, mask) = z_rect_all(i);
end
% Save river bottom elevation for each cross-section
save(fullfile(out_hydro, sprintf('HM_z_rect_all_%s.mat', run_tag)), 'z_rect_all');
% Variables needed by the morphodynamic module (MorphodynamicModule.m)
save(fullfile(out_hydro, sprintf('HM_workspace_%s.mat', run_tag)), ...
    'all_centerlines', 'sections_ids', 'bif_info', 'conf_info', ...
    'widths_all', 'slope_sect', 'z_rect_all', 'h_rect_all', 'Q_map', ...
    'predecessor', 'topo_order', 'network_root', ...
    'x_matrix', 'y_matrix', 'z_matrix');

%% DTM and wet area mask
[mask_wetarea, Rmask_wetarea] = readgeoraster(fullfile(input_dir, wetarea_file));
mask_wetarea = logical(mask_wetarea);
% Read dtm as raster
[raster_dtm, R_dtm] = readgeoraster(fullfile(input_dir, dtm_file));
raster_dtm = double(raster_dtm);
raster_dtm(raster_dtm == -32767) = NaN;
% Are R_dtm and Rmask_wetarea aligned?
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

%% Wet area mask on the DTM grid
% World coordinates of every mask pixel
[cols_mask, rows_mask] = meshgrid(1:size(mask_wetarea,2), 1:size(mask_wetarea,1));
[xw_mask, yw_mask] = Rmask_wetarea.intrinsicToWorld(cols_mask, rows_mask);

% Keep only the TRUE pixels of the mask
mask_true_idx = find(mask_wetarea);
x_mask_true = xw_mask(mask_true_idx);
y_mask_true = yw_mask(mask_true_idx);

% Convert their world coordinates into intrinsic indices of the DTM
[col_in_dtm, row_in_dtm] = R_dtm.worldToIntrinsic(x_mask_true, y_mask_true);

% Round to integer indices
col_in_dtm = round(col_in_dtm);
row_in_dtm = round(row_in_dtm);

% Discard indices outside the DTM
valid = col_in_dtm >= 1 & col_in_dtm <= size(raster_dtm,2) & ...
        row_in_dtm >= 1 & row_in_dtm <= size(raster_dtm,1);

col_in_dtm = col_in_dtm(valid);
row_in_dtm = row_in_dtm(valid);

% Build mask_extended with linear indices
mask_extended = false(size(raster_dtm));
lin_idx = sub2ind(size(raster_dtm), row_in_dtm, col_in_dtm);
mask_extended(lin_idx) = true;
% Segments without discharge (no_flow_segments): their wet pixels are not
% corrected, the raw DTM is kept there
mask_extended = exclude_no_flow_pixels(mask_extended, R_dtm, x_matrix, y_matrix, all_centerlines, no_flow_segments);

% Check: sample some TRUE pixels and verify that the coordinates match
sample_idx = mask_true_idx(1:100:end);  % one TRUE pixel every 100
x_check = xw_mask(sample_idx);
y_check = yw_mask(sample_idx);

[col_check, row_check] = R_dtm.worldToIntrinsic(x_check, y_check);
[x_back, y_back] = R_dtm.intrinsicToWorld(round(col_check), round(row_check));

max_err_x = max(abs(x_check - x_back));
max_err_y = max(abs(y_check - y_back));
fprintf('Max roundtrip error: X=%.4f m, Y=%.4f m\n', max_err_x, max_err_y);

%% Riverbed elevation along cross-sections
% Extract valid points (not NaN) and plot riverbed elevation
mask_valid = ~isnan(z_rect_all_matrix);
x = x_matrix(mask_valid);
y = y_matrix(mask_valid);
z = z_rect_all_matrix(mask_valid);
fig = figure;
scatter(x, y, 10, z, 'filled');
clim_riverbed = prctile(z, [1 99]);   % percentiles instead of min/max: robust to isolated outliers
caxis(clim_riverbed);
ylabel(colorbar, 'Z (m)');
xlabel('X'); ylabel('Y');
title('Riverbed elevation along cross-sections');
axis equal;
save_fig(fig, save_figures, out_fig_hydro, 'HM_riverbed_cross_sections');

%% Interpolation on the DTM grid and correction
% Create Interpolant
F = scatteredInterpolant(x, y, z, 'natural', 'none');

% Create interpolation grid
[cols, rows] = meshgrid(1:size(raster_dtm,2), 1:size(raster_dtm,1));
[xq, yq] = R_dtm.intrinsicToWorld(cols, rows);  % Spatial coordinates
zq = F(xq, yq);  % Interpolation over the grid

% Correction
raster_dtm_corrected = raster_dtm;

% Never raise the bed above the original LiDAR elevation: the model only adds
% bathymetric information (deepening), it never overwrites upwards what the
% LiDAR has already measured.
n_would_raise = nnz(mask_extended & (zq > raster_dtm));

zq_safe = min(zq, raster_dtm);
raster_dtm_corrected(mask_extended) = zq_safe(mask_extended);

fig = figure;
imAlpha = true(size(raster_dtm_corrected));
imAlpha(raster_dtm_corrected == 0 | isnan(raster_dtm_corrected)) = false;
imagesc(R_dtm.XWorldLimits, R_dtm.YWorldLimits, raster_dtm_corrected, 'AlphaData', imAlpha);
axis image;
colormap(parula);
valid_vals = raster_dtm_corrected(~isnan(raster_dtm_corrected));
caxis(prctile(valid_vals, [1 99]));
cb = colorbar;
ylabel(cb, 'Z (m)');
xlabel('X'); ylabel('Y');
title('DTM after bathymetric correction');
save_fig(fig, save_figures, out_fig_hydro, 'HM_DTM_corrected');

%% Bathymetric correction map - DoD
DoD = raster_dtm - raster_dtm_corrected;

fig = figure;
imAlpha = true(size(DoD));
imAlpha(DoD == 0 | isnan(DoD)) = false;   % transparent where no correction was applied

imagesc(R_dtm.XWorldLimits, R_dtm.YWorldLimits, DoD, 'AlphaData', imAlpha);
set(gca, 'Color', 'white');           % white background under the transparent pixels
axis image;
colormap(parula);
nz = DoD(DoD > 0 & ~isnan(DoD));
L = prctile(nz, 98);      % dynamic limit, robust to isolated extreme pixels
caxis([0 L]);
cb = colorbar;
ylabel(cb, 'ΔZ (m)');
title('Bathymetric correction map');
xlabel('X'); ylabel('Y');
save_fig(fig, save_figures, out_fig_hydro, 'HM_DoD');

% Save
geotiffwrite(fullfile(out_hydro, sprintf('HM_DTM_corrected_%s.tif', run_tag)), raster_dtm_corrected, R_dtm, 'CoordRefSysCode', crs_code);
geotiffwrite(fullfile(out_hydro, sprintf('HM_DoD_%s.tif', run_tag)), DoD, R_dtm, 'CoordRefSysCode', crs_code);
disp('Done!');

%% Checking
is_bif_conf = ismember(all_centerlines(:,1), bif_ids) | ismember(all_centerlines(:,1), conf_ids);
z_coords_check = all_centerlines(:,4);
residuo_formula = z_rect_all - (z_coords_check - h_rect_all);
inconsistenti = find(abs(residuo_formula) > 0.01 & ~is_bif_conf);
fprintf('%d points (single channel) where z_rect is NOT consistent with z_coords - h_rect\n', numel(inconsistenti));
if ~isempty(inconsistenti)
    disp(all_centerlines(inconsistenti(1:min(10,end)), [1 6]))
end

% outlier_flags is produced by the block "Clean anomalous single values"
righe_pulite = find(outlier_flags);
righe_inconsistenti = inconsistenti;  % from the previous consistency check

fprintf('Rows cleaned by the anomaly block: %d\n', numel(righe_pulite));
fprintf('Inconsistent rows found now: %d\n', numel(righe_inconsistenti));
fprintf('Overlap: %d\n', numel(intersect(righe_pulite, righe_inconsistenti)));

% The corrected bed elevation must decrease (or stay flat) downstream, with a
% tolerance for local noise: a large, isolated upward jump is suspicious
fprintf('%-8s %-12s %-16s\n', 'Segment', 'N_rises', 'Max_rise_m');
for t = unique(all_centerlines(:,6))'
    idx = find(all_centerlines(:,6) == t);
    if numel(idx) < 3, continue; end
    dz = diff(z_rect_all(idx));
    rises = dz(dz > 1);   % 1 m threshold, adjust if too/not sensitive enough
    if ~isempty(rises)
        fprintf('%-8d %-12d %-16.2f\n', t, numel(rises), max(rises));
    end
end

%% Diagnostics (specific to the 2025 network)
if run_diagnostics
    % Anomaly of segment 12
    idx = find(all_centerlines(:,6) == 12);
    if numel(idx) > 1
        dz = diff(z_rect_all(idx));
        [max_dz, k] = max(dz);
        row_prev = idx(k); row_next = idx(k+1);

        fprintf('Segment 12: max jump = %.2f m between ID %d and ID %d\n', max_dz, ...
            all_centerlines(row_prev,1), all_centerlines(row_next,1));

        % Context: a few sections before and after the jump
        ctx = max(1,k-3):min(numel(idx),k+3);
        T = table(all_centerlines(idx(ctx),1), widths_all(idx(ctx)), sect_slope(idx(ctx)), ...
            Q_map(idx(ctx)), h_rect_all(idx(ctx)), z_rect_all(idx(ctx)), ...
            ismember(all_centerlines(idx(ctx),1), bif_ids), ...
            ismember(all_centerlines(idx(ctx),1), conf_ids), ...
            'VariableNames', {'ID','width','slope','Q','h_rect','z_rect','is_bif','is_conf'});
        disp(T)
    end

    % Do the segments oscillate? Raw z_DTM (z_original, before the glitch
    % correction) of the segments that oscillated before the nodata fix
    % (19/11/2026) and of those with the largest raw jumps.
    %   roughness  standard deviation of z around its moving median on 7 sections
    %   Jump_raw   largest jump between consecutive sections (a step is not a glitch)
    %   glitches   sections corrected by the glitch rule, largest |excess|
    % Outcome: 'no oscillation' if roughness < rough_tol, otherwise
    % 'TO BE INSPECTED MANUALLY' (an oscillation over several sections is not
    % removed by the glitch rule).
    rough_tol = 0.20;   % [m]
    tronchi_sospetti = [0 12 19 26 30 32 34 35 37 38 39];
    fprintf('\nSegments checked for oscillations of z_DTM (raw z)\n');
    fprintf('%-8s %-11s %-10s %-9s %-12s %s\n', 'Segment', 'roughness', 'Jump_raw', 'glitches', 'max|excess|', 'Outcome');
    for t = tronchi_sospetti
        idx_t = find(all_centerlines(:,6) == t);
        if numel(idx_t) < 7, continue; end
        z0 = z_original(idx_t);
        rough = std(z0 - movmedian(z0, 7), 'omitnan');
        ng = nnz(outlier_flags(idx_t));
        mx = max([0; abs(glitch_excess(idx_t(outlier_flags(idx_t))))]);
        if rough < rough_tol
            esito = 'no oscillation';
        else
            esito = 'TO BE INSPECTED MANUALLY (oscillation)';
        end
        fprintf('%-8d %-11.3f %-10.2f %-9d %-12.2f %s\n', t, rough, max(abs(diff(z0))), ng, mx, esito);
    end

    for t = [32 35]
        idx = find(all_centerlines(:,6) == t);
        if numel(idx) < 2, continue; end
        dz = diff(z_original(idx));   % raw z_coords, before the glitch correction
        [max_dz, k] = max(abs(dz));
        row_prev = idx(k); row_next = idx(k+1);
        fprintf('\nSegment %d: max jump = %.2f m between ID %d and ID %d\n', t, dz(k), ...
            all_centerlines(row_prev,1), all_centerlines(row_next,1));
        ctx = max(1,k-3):min(numel(idx),k+3);
        T = table(all_centerlines(idx(ctx),1), widths_all(idx(ctx)), sect_slope(idx(ctx)), ...
            Q_map(idx(ctx)), h_rect_all(idx(ctx)), z_rect_all(idx(ctx)), ...
            ismember(all_centerlines(idx(ctx),1), bif_ids), ...
            ismember(all_centerlines(idx(ctx),1), conf_ids), ...
            'VariableNames', {'ID','width','slope','Q','h_rect','z_rect','is_bif','is_conf'});
        disp(T)
    end

    % ID 168 is the start of segment 8: check whether it is registered as a confluence node
    match = find(conf_ids == 168);
    if ~isempty(match)
        branch_ids = conf_info{match,2};
        fprintf('Confluence at node %d: upstream branches = %s\n', 168, mat2str(branch_ids(:)'));
        for b = branch_ids(:)'
            tronco_b = all_centerlines(all_centerlines(:,1)==b, 6);
            fprintf('  Branch ID %d -> segment %d\n', b, tronco_b);
        end
    else
        fprintf('ID 168 not found in conf_info as an explicit confluence node.\n');
    end

    % Pending checks
    for t = [6 7]
        row_start = find(all_centerlines(:,6)==t, 1, 'first');
        if isempty(row_start), continue; end
        seg_idx = find(sections_ids(:,1) == all_centerlines(row_start,1));
        up = predecessor{seg_idx};
        if isempty(up)
            fprintf('Segment %d: NO predecessor (actual network head)\n', t);
        else
            tronco_predecessore = codes_by_section(up);
            fprintf('Segment %d: predecessor = row %d -> SEGMENT CODE %d\n', t, up, tronco_predecessore);
        end
    end
    % Follow the chain of predecessors up to the actual network head, printing Q at each step
    t = 7;
    while any(all_centerlines(:,6)==t)
        row_start = find(all_centerlines(:,6)==t, 1, 'first');
        Q_t = Q_map(row_start);
        seg_idx = find(sections_ids(:,1) == all_centerlines(row_start,1));
        up = predecessor{seg_idx};
        if isempty(up)
            fprintf('Segment %d (NETWORK HEAD): Q = %.4f\n', t, Q_t);
            break
        end
        t_up = codes_by_section(up(1));
        fprintf('Segment %d: Q = %.4f  <-  predecessor segment %d\n', t, Q_t, t_up);
        t = t_up;
    end
end

%% ============================================================
%% VALIDATION
%% ============================================================
field_pts = read_field_points(fullfile(input_dir, gnss_file), fullfile(input_dir, adcp_file), ...
    gnss_station_split, adcp_min_spacing);   % GNSS cross-sections + ADCP transects
[Tstats, Tsec] = validation_gnss_mean(field_pts, raster_dtm, raster_dtm_corrected, R_dtm, mask_extended, out_val)
% RMSE/MAE improve (the model corrects the large, systematic errors well),
% but it introduces its own "typical" scatter, which the rectangularisation alone does not explain.

if run_diagnostics
    x0 = 406750; y0 = 5010300;   % approximate coordinates of station s006
    dmat = sqrt((x_matrix - x0).^2 + (y_matrix - y0).^2);
    hit_rows = find(min(dmat, [], 2, 'omitnan') < 60);
    T = table(all_centerlines(hit_rows,1), all_centerlines(hit_rows,6), all_centerlines(hit_rows,4), ...
        widths_all(hit_rows), sect_slope(hit_rows), Q_map(hit_rows), h_rect_all(hit_rows), z_rect_all(hit_rows), ...
        'VariableNames', {'ID','segment','z_coords','width','slope','Q','h_rect','z_rect'});
    disp(T)
end

disp('Preprocessing & Hydraulic module: all outputs saved in Output_folder.');

%% Local functions
function save_fig(fig, do_save, folder, name)
% Save a figure as PNG in folder (if do_save) and close it when running without display
if do_save
    exportgraphics(fig, fullfile(folder, [name '.png']), 'Resolution', 150);
end
if ~usejava('desktop')
    close(fig);
end
end
