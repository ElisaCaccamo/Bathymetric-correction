function [W_bf, left_pt, right_pt, A_bf] = W_bankfull_from_dtm( ...
    i, z_rect_all, h_bf, ...
    Z, Xgrid, Ygrid, all_centerlines, opts)

% INPUTS:
% i                  : current section index (== all_centerlines(i,1), IDs are dense 1..N)
% z_rect_all         : vector [n x 1] riverbed elevation
% h_bf               : bankfull water depth for section i [m]
% Xgrid, Ygrid, Z    : DTM grid compatible with interp2
% all_centerlines    : [N x 5] id, x, y, z (free surface), width_old
% opts               : struct with optional fields
%      - ds              : sampling step along the transect (m)
%      - max_half_len    : (default 100) max half-length of the transect (m)
%      - interp_method    : (default 'linear') interp2 method
%      - h_bf_all         : (optional vector) bankfull depths already computed,
%                           used as fallback for the previous section
%      - verbose          : (default false) print diagnostics
%
% OUTPUTS:
%  W_bf [1x1], A_bf [1x1], left_pt (x,y) [1x2], right_pt (x,y) [1x2]

% --- default outputs (used only if every path below fails) ---
W_bf = NaN;
left_pt = [NaN, NaN];
right_pt = [NaN, NaN];
A_bf = NaN;

x_coords = all_centerlines(:, 2);
y_coords = all_centerlines(:, 3);
W_old = all_centerlines(:, 5);
N = numel(x_coords);
ds = opts.ds;
max_half = opts.max_half_len;
interp_method = opts.interp_method;
verbose = isfield(opts, 'verbose') && opts.verbose;

if isfield(opts,'h_bf_all') && ~isempty(opts.h_bf_all)
    h_bf_all = opts.h_bf_all;
else
    h_bf_all = [];
end

s = -max_half:ds:max_half; % transect positions relative to center

xc = x_coords(i); yc = y_coords(i);
ht = h_bf;
z_target = z_rect_all(i) + ht;

% --- tangent / normal to the centerline ---
if i == 1
    dx_t = x_coords(i+1) - xc;
    dy_t = y_coords(i+1) - yc;
elseif i == N
    dx_t = xc - x_coords(i-1);
    dy_t = yc - y_coords(i-1);
else
    dx_t = x_coords(i+1) - x_coords(i-1);
    dy_t = y_coords(i+1) - y_coords(i-1);
end
tvec = [dx_t, dy_t];
if norm(tvec) == 0
    warning('Section %d: undefined centerline direction (coincident points).', i);
    return; % direzione non definita
end
tvec = tvec / norm(tvec);     % tangente normalizzata
nvec = [-tvec(2), tvec(1)];   % normale (perpendicolare)

% --- transect sampling ---
xsec = xc + nvec(1) * s;
ysec = yc + nvec(2) * s;
zsec = interp2(Xgrid, Ygrid, Z, xsec, ysec, interp_method, NaN);

if verbose
    fprintf('Section %d diagnostics:\n', i);
    fprintf('  s range = [%.2f, %.2f], max_half=%.2f, ds=%.3f\n', min(s), max(s), max_half, ds);
    fprintf('  xsec range = [%.3f, %.3f], Xgrid range = [%.3f, %.3f]\n', min(xsec), max(xsec), min(Xgrid(:)), max(Xgrid(:)));
    fprintf('  ysec range = [%.3f, %.3f], Ygrid range = [%.3f, %.3f]\n', min(ysec), max(ysec), min(Ygrid(:)), max(Ygrid(:)));
    nz = zsec(~isnan(zsec));
    if isempty(nz)
        fprintf('  zsec: ALL NaN -> interp2 out of grid or CRS mismatch\n');
    else
        fprintf('  zsec min=%.3f max=%.3f nNaN=%d\n', min(nz), max(nz), sum(isnan(zsec)));
    end
    fprintf('  z_rect_all(i)=%.3f, ht=%.3f, z_target=%.3f\n', z_rect_all(i), ht, z_target);

    in_grid_mask = xsec >= min(Xgrid(:)) & xsec <= max(Xgrid(:)) & ...
                   ysec >= min(Ygrid(:)) & ysec <= max(Ygrid(:));
    n_out_of_grid = sum(~in_grid_mask);
    fprintf('  Section %d: %d / %d points out of grid\n', i, n_out_of_grid, numel(xsec));
end
if any(isnan(xsec)) || any(isnan(ysec))
    warning('Section %d: xsec or ysec contain NaN values.', i);
end

fallback_code = zeros(N,1);
% 0 = normal crossing, 1 = exact-zero sample, 2 = clamp used,
% 3 = percentile depression used, 4 = W_old used,
% 5 = halfW relaxed after clamp, 6 = mirrored single-side

% --- crossings between zsec and z_target (fallback A-C) ---
D = zsec - z_target;
valid = ~isnan(D);
idx_valid = find(valid);
if numel(idx_valid) < 2
    warning('Section %d: insufficient valid samples on transect (nValid=%d).', i, numel(idx_valid));
    return;
end
s_valid = s(idx_valid);
D_valid = D(idx_valid);
zsec_valid = zsec(idx_valid);

tol = 1e-3;        % tolerance for "near-zero" crossing [m]
eps_clamp = 1e-2;  % clamp margin when z_target outside profile [m]  (kept for reference, not otherwise used)
p_quant = 0.10;    % percentile for depression heuristic (10%)

zmin = min(zsec_valid);
zmax = max(zsec_valid);

% 1) standard zero-crossings with tolerance
pairs = find( D_valid(1:end-1).*D_valid(2:end) < 0 | abs(D_valid(1:end-1))<tol | abs(D_valid(2:end))<tol );
cross_positions = [];
for k = pairs
    s1 = s_valid(k); s2 = s_valid(k+1);
    d1 = D_valid(k); d2 = D_valid(k+1);
    if abs(d2 - d1) < eps
        s_cross = 0.5*(s1 + s2);
    else
        s_cross = s1 + (0 - d1) * (s2 - s1) / (d2 - d1);
    end
    cross_positions(end+1) = s_cross; %#ok<AGROW>
end
zero_idx = find(abs(D_valid) < tol);
if ~isempty(zero_idx)
    cross_positions = unique([cross_positions, s_valid(zero_idx)']);
    fallback_code(i) = 1;
end

% 2) if no crossings, clamp z_target to profile extremes and retry
clamp_attempted = false;
if isempty(cross_positions)
    if z_target < zmin
        z_target_try = zmin;
        clamp_attempted = true;
    elseif z_target > zmax
        z_target_try = zmax;
        clamp_attempted = true;
    else
        z_target_try = z_target;
    end

    if clamp_attempted
        D_try = zsec_valid - z_target_try;
        pairs_try = find( D_try(1:end-1).*D_try(2:end) < 0 | abs(D_try(1:end-1))<tol | abs(D_try(2:end))<tol );
        for k = pairs_try
            s1 = s_valid(k); s2 = s_valid(k+1);
            d1 = D_try(k); d2 = D_try(k+1);
            if abs(d2 - d1) < eps
                s_cross = 0.5*(s1 + s2);
            else
                s_cross = s1 + (0 - d1) * (s2 - s1) / (d2 - d1);
            end
            cross_positions(end+1) = s_cross; %#ok<AGROW>
        end
        zero_idx_try = find(abs(D_try) < tol);
        if ~isempty(zero_idx_try)
            cross_positions = unique([cross_positions, s_valid(zero_idx_try)']);
        end
        if ~isempty(cross_positions)
            z_target = z_target_try; % accept clamped target going forward
            fallback_code(i) = 2;
            warning('Section %d: z_target clamped to within profile range (clamp used).', i);
        else
            warning('Section %d: clamp of z_target attempted but no crossings found.', i);
        end
    end
end

% 3) if still no crossings, detect channel depression via percentile/min
s_left = NaN; s_right = NaN;
if isempty(cross_positions)
    center_hold = min(max_half, 0.75 * max_half);
    center_mask = abs(s_valid) <= center_hold;
    if any(center_mask)
        z_center = zsec_valid(center_mask);
        s_center = s_valid(center_mask); %#ok<NASGU>
        zq = quantile(z_center, p_quant);
        low_idx_local = find(center_mask & (zsec_valid <= zq));
        if isempty(low_idx_local)
            [~, imin] = min(abs(s_valid));
            window_idx = max(1, imin-3):min(numel(s_valid), imin+3);
            [~, minpos] = min(zsec_valid(window_idx));
            idx_global = window_idx(minpos);
            if s_valid(idx_global) < 0, s_left = s_valid(idx_global); else, s_right = s_valid(idx_global); end
        else
            s_low = s_valid(low_idx_local);
            left_candidates  = s_low(s_low < 0);
            right_candidates = s_low(s_low > 0);
            if ~isempty(left_candidates)
                s_left = max(left_candidates);
            end
            if ~isempty(right_candidates)
                s_right = min(right_candidates);
            end
        end
        if ~isnan(s_left) || ~isnan(s_right)
            fallback_code(i) = 3;
            warning('Section %d: banks inferred from low-percentile depression (percentile=%g).', i, p_quant);
        end
    end
end

% 4) if still no crossings, use W_old to build banks
if isempty(cross_positions)
    if ~isnan(W_old(i))
        W_bf = W_old(i);
        s_left = -0.5 * W_bf;
        s_right = 0.5 * W_bf;
        A_bf = W_bf * ht; % approximate area
        fallback_code(i) = 4;
        warning('Section %d: fallback D used (W_old = %.3f m).', i, W_bf);
        left_pt = [xc + nvec(1)*s_left, yc + nvec(2)*s_left];
        right_pt = [xc + nvec(1)*s_right, yc + nvec(2)*s_right];
        return;
    end
end

% 5) if standard crossing(s) were found, apply halfW filter
if ~isempty(cross_positions)
    halfW = NaN;
    if ~isnan(W_old(i))
        halfW = 0.5 * W_old(i);
    end

    if isnan(halfW)
        filtered_cross = cross_positions;
    else
        filtered_cross = cross_positions(abs(cross_positions) >= halfW);
    end

    if isempty(filtered_cross) && clamp_attempted
        filtered_cross = cross_positions; % ignore halfW after clamp
        fallback_code(i) = 5;
        warning('Section %d: halfW filter relaxed/ignored because clamp was attempted.', i);
    end

    if ~isempty(filtered_cross)
        left_candidates  = filtered_cross(filtered_cross < 0);
        right_candidates = filtered_cross(filtered_cross > 0);
        if ~isempty(left_candidates)
            s_left = max(left_candidates);
        end
        if ~isempty(right_candidates)
            s_right = min(right_candidates);
        end
        if fallback_code(i) == 0 && ~isempty(zero_idx)
            fallback_code(i) = 1;
        end
    else
        if ~isempty(cross_positions)
            warning('Section %d: no crossing satisfies abs(s) >= halfW (halfW=%.2f) and filter not relaxed.', i, halfW);
        end
    end
end

% 6) if one bank still undefined, mirror the other
if (isnan(s_left) && ~isnan(s_right)) || (~isnan(s_left) && isnan(s_right))
    if isnan(s_left) && ~isnan(s_right)
        s_left = -abs(s_right);
    elseif ~isnan(s_left) && isnan(s_right)
        s_right =  abs(s_left);
    end
    if abs(s_left) > max_half
        s_left = -max_half;
    end
    if abs(s_right) > max_half
        s_right = max_half;
    end
    fallback_code(i) = 6;
    warning('Section %d: single-side detected -> mirrored found side (s_left=%.3f, s_right=%.3f).', i, s_left, s_right);
end

if isnan(s_left) || isnan(s_right)
    warning('Section %d: could not identify both banks after fallbacks (s_left=%g, s_right=%g). Skipping section.', i, s_left, s_right);
    return;
end

% --- bank coordinates ---
xL = xc + nvec(1) * s_left;  yL = yc + nvec(2) * s_left;
xR = xc + nvec(1) * s_right; yR = yc + nvec(2) * s_right;

% --- cross-sectional area (integrate z_target - zsec between banks) ---
nSamples = max(200, round((s_right - s_left) / mean(diff(s_valid))));
s_fine = linspace(s_left, s_right, nSamples);
zsec_fine = interp1(s_valid, zsec_valid, s_fine, 'linear', NaN);

if any(isnan(zsec_fine))
    warning('Section %d: NaN values while interpolating the cross-section profile, skip.', i);
    return;
end

area_section = trapz(s_fine, (z_target - zsec_fine)); % m^2

if area_section <= 0
    % --- fallback: try previous section's z_target/h_bf ---
    warning('Section %d: non-positive area (%.3f), trying fallback with previous z_target.', i, area_section);
    if verbose
        fprintf('ht (h_bf) = %.3f m, current z_target = %.3f\n', ht, z_target);
    end
    if i > 1
        z_target_prev = NaN;
        if ~isempty(h_bf_all) && numel(h_bf_all) >= i-1 && ~isnan(h_bf_all(i-1))
            ht_prev = h_bf_all(i-1);
            z_target_prev = z_rect_all(i-1) + ht_prev;
            info_str = sprintf('ht_prev = %.3f (global h_bf_all)', ht_prev);
        else
            z_target_prev = z_rect_all(i-1) + ht;
            info_str = sprintf('using current ht = %.3f with previous bed elevation', ht);
        end

        area_prev = trapz(s_fine, (z_target_prev - zsec_fine));
        if verbose
            fprintf('z_target_prev = %.3f (%s), area_prev = %.3f\n', z_target_prev, info_str, area_prev);
        end

        if area_prev > 0
            area_section = area_prev;
            W_from_area = area_section / ht;
            if ~isnan(W_old(i)) && (W_from_area < W_old(i))
                W_bf = NaN;
                warning('Section %d: W_from_area fallback (%.2f m) < W_old (%.2f m) -> NaN', i, W_from_area, W_old(i));
                return;
            else
                W_bf = W_from_area;
                left_pt = [xL, yL];
                right_pt = [xR, yR];
                A_bf = area_section;
                return;
            end
        end
        % previous-area fallback also failed -> use W_old
        if ~isnan(W_old(i))
            W_bf = W_old(i);
            left_pt = [xc - 0.5*W_bf * nvec(1), yc - 0.5*W_bf * nvec(2)];
            right_pt = [xc + 0.5*W_bf * nvec(1), yc + 0.5*W_bf * nvec(2)];
            A_bf = W_bf * ht;
            warning('Section %d: fallback with W_old -> W_bf = %.3f (estimated A_bf = %.3f).', i, W_bf, A_bf);
            return;
        else
            warning('Section %d: previous-area fallback failed and W_old missing -> skip', i);
            return;
        end
    else
        warning('Section %d: non-positive area at i==1, no previous section available -> skip', i);
        return;
    end
end

% --- NORMAL PATH (area_section > 0): this used to be unreachable due to a
%     missing 'end' in the original file. Fixed here. ---
W_from_area = area_section / ht;
W_bf = W_from_area;
left_pt = [xL, yL];
right_pt = [xR, yR];
A_bf = area_section;

% Optional conservative rule (uncomment to enforce W_bf >= W_old):
% if ~isnan(W_old(i)) && (W_from_area < W_old(i))
%     W_bf = NaN;
%     warning('Section %d: W_from_area (%.2f m) < W_old (%.2f m) -> NaN', i, W_from_area, W_old(i));
% end

end
