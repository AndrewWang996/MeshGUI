function logfz = getLogFz(meshname, whichKeyframe)
    keyframe = getKeyframe(meshname, whichKeyframe);
    fz = keyframe.fz;
    
    [endnodes, weights, predecessors] = getSpanningTree(meshname);
    tree = graph(endnodes(:,1), endnodes(:,2), weights);
    
    anchorIndices = getAnchorIndices(meshname);
    anchorIndex = anchorIndices(1);
    
    
%     g = complex( log(abs(allFz)), angle(allFz(end)) + cumsum(angle(allFz./circshift(allFz, 1))) );
    predecessors(anchorIndex) = anchorIndex;
    individualAngles = angle( fz ./ fz(predecessors) );
    accumulatedAngles = accumulateAlongEdges(tree, anchorIndex, individualAngles);
    logfz = complex( log(abs(fz)), angle(fz(anchorIndex)) + accumulatedAngles );
    
end