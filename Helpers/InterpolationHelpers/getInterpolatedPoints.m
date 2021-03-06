function points = getInterpolatedPoints(nVertices, nKeyframes, coeffs, numTimesPerInterval)
%% A horribly written piece of code that interpolates using Hermite splines
% with input values and derivatives
% numTimesPerInterval refers to how many time points we want to sample from
% each interval.
    nKeyframePairs = nKeyframes - 1;
    
    times = linspace(0, 1, 1 + numTimesPerInterval);
    times = times(1:end-1);
    
    timeMatrix = zeros(4, numTimesPerInterval);
    timeMatrix(4,:) = ones(1, numTimesPerInterval);
    for i=3:-1:1
        timeMatrix(i,:) = timeMatrix(i+1,:) .* times;
    end
    
    points_stack = coeffs * timeMatrix;
    points = zeros(nVertices, nKeyframePairs * numTimesPerInterval + 1);
    for interval = 1 : nKeyframePairs
        startTime = (interval - 1) * numTimesPerInterval + 1;
        endTime = interval * numTimesPerInterval;
        
        startVert = (interval - 1) * nVertices + 1;
        endVert = interval * nVertices;
        points(:, startTime:endTime ) = points_stack(startVert:endVert, :);
    end
    points(:,end) = sum( coeffs(end-nVertices+1:end,:), 2);
end


