function control_points = getControlPoints(interpolated_points)
    control_points = zeros( size(interpolated_points) );
    control_points(1,:) = interpolated_points(1,:);
    control_points(end,:) = interpolated_points(end,:);

    
    num_points = size(interpolated_points, 1);
    tmp_ones = ones( 1, num_points - 3 );
    diag_matrix = 4 * eye(num_points - 2) + diag(tmp_ones, 1) + diag(tmp_ones, -1);
    
    vector = 6 * interpolated_points(2:end-1,:);
    vector(1,:) = vector(1,:) - interpolated_points(1,:);
    vector(end,:) = vector(end,:) - interpolated_points(end,:);
    
    intermediate_control_points = diag_matrix \ vector;
    control_points(2:end-1,:) = intermediate_control_points(:,:);
end