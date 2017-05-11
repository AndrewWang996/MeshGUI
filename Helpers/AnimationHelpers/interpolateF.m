function interpF = interpolateF(info, numTimesPerInterval)
    anchorIndex = info.anchorIndex;
    vertices = info.vertices;
    allVertices = info.allVertices;
    logFz = info.allLogFz;
    allEta = info.allEta;
    all_deta_dt = info.all_deta_dt;
    endNodes = info.endNodes;
    
    tree = graph(endNodes(:,1), endNodes(:,2));
    numKeyframes = size(allEta, 2);
    
    complexVertices = complex(vertices(:,1), vertices(:,2));
    edgeVectors = complexVertices(endNodes(:,2)) - complexVertices(endNodes(:,1));

    % 2) interpolate fz, eta, fzbar
    interpFz = exp( getInterpolatedPointsBezier( logFz, numTimesPerInterval ) );
    interpEta = getInterpolatedPointsHermite(allEta, all_deta_dt, numTimesPerInterval);
    interpFzBar = interpEta ./ conj(interpFz);

    % 3) integrate fz -> Phi, fzbar -> Psi by collecting edge
    % differences.
    edgeDifferencesFz = (edgeVectors / 2) .* (interpFz(endNodes(:,1),:) + interpFz(endNodes(:,2),:));
    edgeDifferencesFzBar = (edgeVectors / 2) .* (interpFzBar(endNodes(:,1),:) + interpFzBar(endNodes(:,2),:));


    % Traverse the graph, accumulating edge values in Phi, Psi
    allTimes = linspace(0, 1, 1 + (numKeyframes - 1) * numTimesPerInterval);
    weight = linearWeight( numKeyframes, allTimes );
    mixedF = allVertices * weight;
    Phi = accumulateAlongEdges(...
        tree,...
        anchorIndex,...
        edgeDifferencesFz,...
        mixedF(anchorIndex,:)... % Set Phi, Psi anchor values as defined in BDHI
    );
    PsiBar = accumulateAlongEdges(...
        tree,...
        anchorIndex,...
        edgeDifferencesFzBar...
    );

    % 4) sum them together
    interpF = Phi + PsiBar;
end

