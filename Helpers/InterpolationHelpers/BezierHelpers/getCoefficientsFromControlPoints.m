function coeffs = getCoefficientsFromControlPoints(control_points)
    S = zeros( size(control_points) );
    S(:,1) = control_points(:,1);
    S(:,end) = control_points(:,end);
    
    S(:,2:end-1) = (1/6) * (...
        control_points(:,1:end-2) ...
        + 4 * control_points(:,2:end-1) ...
        + control_points(:,3:end) ...
    );

    B = control_points;
    
    
    nVertices = size(control_points, 1);
    nKeyframes = size(control_points, 2);
    nSplines = nKeyframes - 1;
    coeffs = zeros( 4 , nVertices * nSplines );
    
    for i = 1:nSplines
        control_matrix = horzcat( ...
            S(:,i), ...
            (1/3) * (2 * B(:,i) +     B(:,i+1)), ...
            (1/3) * (    B(:,i) + 2 * B(:,i+1)), ...
            S(:,i+1) ...
        );
        block = getBezierCoefficients( transpose( control_matrix )  );
        coeffs(1:4, nVertices*(i-1)+1 : nVertices*i) = block;
    end
    
    coeffs = transpose(coeffs);
    
end


function sub_coeffs = getBezierCoefficients(control_points)
% control_points are an (4 x n) complex matrix
% we output an (4 x n) matrix of coefficients (a_i, b_i, c_i, d_i)
% such that P(i,t) = a_i * t^3 + b_i * t^2 + c_i * t + d_i
    
    p0 = control_points(1,:);
    p1 = control_points(2,:);
    p2 = control_points(3,:);
    p3 = control_points(4,:);
    
    a = p3 - 3 * p2 + 3 * p1 -     p0;
    b =      3 * p2 - 6 * p1 + 3 * p0;
    c =               3 * p1 - 3 * p0;
    d =                            p0;
    
    sub_coeffs = vertcat(a,b,c,d);
end

