function points = getInterpolatedPointsBezier(interpolated_points, numTimesPerInterval)
    nVertices = size(interpolated_points, 1);
    nKeyframes = size(interpolated_points, 2);
    
    
    control_points = getControlPoints(interpolated_points);
    
    coeffs = getBezierCoefficients(control_points);

end

