function getInterpolatedPointsTest()
    % test no basic errors
    %gIP_mesh();
    gIP_basic();
end


function gIP_basic()
    allEta = [1 2 3];
    allDeta_dt = [1 0 8];
    numFrames = length(allEta);
    numTimesPerInterval = 100;
    points = getInterpolatedPoints(allEta, allDeta_dt, numTimesPerInterval);
    plot(linspace(0, numFrames - 1,numTimesPerInterval * (numFrames-1) + 1), points);
%     disp(points);
end

function gIP_mesh()
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