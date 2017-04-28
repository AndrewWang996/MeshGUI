function getInterpolatedPoints()
    meshname = 'vert_bar';
    [V,F] = getMesh(meshname);
    
    numKeyframes = countKeyframes(meshname);
    
    allFz = zeros(size(V,1), numKeyframes);
    allFzbar = zeros(size(V,1), numKeyframes);
    
    for whichKeyframe = 1:numKeyframes
        keyframe = getKeyframe(meshname, whichKeyframe);
        allFz(:,whichKeyframe) = keyframe.fz;
        allFzbar(:,whichKeyframe) = keyframe.fzbar;
    end
    
    allEta = conj(allFz) .* allFzbar;
    
    allDeta_dt = zeros( size(allEta) );
    
    coeffs = getHermiteCoefficients(allEta, allDeta_dt);
    points = getInterpolatedPoints(allEta, allDeta_dt, 10);
end