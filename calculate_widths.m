function [widths, sect_x, sect_y, sect_z, z_sect_mean] = calculate_widths(x_interp, y_interp, z_interp, x_world, y_world, river_mask, dtm)
    widths = NaN(size(x_interp));
    step = 0.5; % Grid resolution
    max_extent = 200; % Maximum extent [m] of the orthogonal vector
    t = -max_extent:step:max_extent;
    ortho_vector = NaN(length(x_interp), 2); 

    % Matrixes for cross-sections coordinates 
    sect_x = NaN(length(x_interp), length(t));
    sect_y = NaN(length(x_interp), length(t));
    sect_z = NaN(length(x_interp), length(t)); 

    % Vector for mean elevations per cross-section
    z_sect_mean = NaN(size(x_interp));

    for i = 1:length(x_interp)
        % Centerline direction
        if i == 1
            dx = x_interp(i+1) - x_interp(i);
            dy = y_interp(i+1) - y_interp(i);
        elseif i == length(x_interp)
            dx = x_interp(i) - x_interp(i-1);
            dy = y_interp(i) - y_interp(i-1);
        else
            dx = x_interp(i+1) - x_interp(i-1);
            dy = y_interp(i+1) - y_interp(i-1);
        end

        direction = [dx, dy] / norm([dx, dy]);
        ortho = [-direction(2), direction(1)];
        ortho_vector(i, :) = ortho;

        x_cross = x_interp(i) + t * ortho(1);
        y_cross = y_interp(i) + t * ortho(2);

        % Check with limits
        valid_points = x_cross >= min(x_world) & x_cross <= max(x_world) & ...
                       y_cross >= min(y_world) & y_cross <= max(y_world);

        x_cross = x_cross(valid_points);
        y_cross = y_cross(valid_points);

        % In-mask check
        in_mask = interp2(x_world, y_world, river_mask, x_cross, y_cross, 'nearest', 0);
        if isempty(in_mask) || all(in_mask == 0)
            continue;
        end

        idx_start = find(in_mask > 0, 1, 'first');
        idx_end = find(in_mask > 0, 1, 'last');
        if isempty(idx_start) || isempty(idx_end)
            continue;
        end

        % Calculate width
        widths(i) = sqrt((x_cross(idx_end) - x_cross(idx_start))^2 + ...
                         (y_cross(idx_end) - y_cross(idx_start))^2);

        dim = abs(idx_end - idx_start) + 1;

        % Coordinates
        x_section = x_cross(idx_start:idx_end);
        y_section = y_cross(idx_start:idx_end);

        sect_x(i, 1:dim) = x_section;
        sect_y(i, 1:dim) = y_section;

        % Interpolate elevation for each cross section point from DTM
        z_cross = interp2(x_world, y_world, dtm, x_cross(idx_start:idx_end), y_cross(idx_start:idx_end), 'linear');
        
        % Take out negative values from mean calculation
        z_valid = z_cross(z_cross >= 0);
        
        if ~isempty(z_valid)
            z_sect_mean(i) = mean(z_valid);
        else
            z_sect_mean(i) = z_interp(i);  
        end
        
        % Assign mean elevation to all points along the cross-section
        sect_z(i, 1:dim) = z_sect_mean(i);
    end

    % Edge case handling
    for i = 1:length(x_interp)
        if i == 1
            widths(i) = widths(i+1); 
            z_sect_mean(i) = z_sect_mean(i+1);
        elseif i == length(x_interp)
            widths(i) = widths(i-1); 
            z_sect_mean(i) = z_sect_mean(i-1);
        end

        ortho = ortho_vector(i, :);
        L = widths(i) / 2;
        t_local = -L:step:L;

        x_cross = x_interp(i) + t_local * ortho(1);
        y_cross = y_interp(i) + t_local * ortho(2);

        sect_x(i, 1:length(t_local)) = x_cross;
        sect_y(i, 1:length(t_local)) = y_cross;

        % Assign mean elevation to the edges
        sect_z(i, 1:length(t_local)) = z_sect_mean(i);
    end
    z_sect_mean = z_sect_mean(:); % vettore colonna
end