function points = getInterpolatedPointsHermite(control_points, tangents, numTimesPerInterval)
    nVertices = size(control_points, 1);
    nKeyframes = size(control_points, 2);
    
    coeffs = getHermiteCoefficients(control_points, tangents);
    
    points = getInterpolatedPoints(nVertices, nKeyframes, coeffs, numTimesPerInterval);
end

