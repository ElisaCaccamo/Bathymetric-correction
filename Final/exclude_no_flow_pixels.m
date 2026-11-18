function [mask_out, n_excluded] = exclude_no_flow_pixels(mask_in, R, x_matrix, y_matrix, all_centerlines, no_flow_segments)
% EXCLUDE_NO_FLOW_PIXELS  Remove from the correction mask the wet pixels that
% belong to segments without discharge at the flight (no_flow_segments).
%
% The bed of those segments is not computed (Q = 0, no uniform-flow depth),
% so their sections are left out of the bed interpolation. The interpolation
% of the neighbouring sections would still reach their pixels, so the pixels
% are removed from the mask and keep the raw DTM.
%
% Each wet pixel is assigned to the nearest cross-section point (any section
% of any segment); a pixel is removed when that point belongs to a segment in
% no_flow_segments.
%
% INPUTS
%   mask_in           logical mask of the corrected area (DTM grid)
%   R                 spatial reference of the DTM grid
%   x_matrix, y_matrix  coordinates of the cross-section points, one row per section
%   all_centerlines   [ID X Y z_DTM W segment_code], one row per section
%   no_flow_segments  segment codes without discharge (empty: nothing is removed)
%
% OUTPUTS
%   mask_out          mask without the pixels of the segments without discharge
%   n_excluded        number of removed pixels

mask_out = mask_in;
n_excluded = 0;
if isempty(no_flow_segments) || ~any(ismember(all_centerlines(:,6), no_flow_segments))
    return
end

% Cross-section points and the code of their segment
seg_code = repmat(all_centerlines(:,6), 1, size(x_matrix, 2));
ok = isfinite(x_matrix) & isfinite(y_matrix);
is_no_flow = double(ismember(seg_code(ok), no_flow_segments));

% Nearest cross-section point of every wet pixel (duplicated points are merged
% by scatteredInterpolant, which then averages the flag: > 0.5 = no flow)
w_state = warning('off', 'MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
F = scatteredInterpolant(x_matrix(ok), y_matrix(ok), is_no_flow, 'nearest', 'nearest');
warning(w_state);
[r, c] = find(mask_in);
[xw, yw] = intrinsicToWorld(R, c, r);
excl = F(xw, yw) > 0.5;

mask_out(sub2ind(size(mask_in), r(excl), c(excl))) = false;
n_excluded = nnz(excl);
fprintf('Segments without discharge (%s): %d wet pixels left uncorrected (%.1f%% of the wet area)\n', ...
    num2str(no_flow_segments(:)'), n_excluded, 100 * n_excluded / max(numel(r), 1));
end
