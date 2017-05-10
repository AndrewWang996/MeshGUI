function coeffs = getBezierCoefficients(control_points)
% control_points are an (n x 4) complex matrix
% we output an (n x 4) matrix of coefficients (a_i, b_i, c_i, d_i)
% such that P(i,t) = a_i * t^3 + b_i * t^2 + c_i * t + d_i
    
    p0 = control_points(1,:);
    p1 = control_points(2,:);
    p2 = control_points(3,:);
    p3 = control_points(4,:);
    
    a = p3 - 3 * p2 + 3 * p1 -     p0;
    b =      3 * p2 - 6 * p1 + 3 * p0;
    c =               3 * p1 - 3 * p0;
    d =                            p0;
    
    coeffs = horzcat(a,b,c,d);
end