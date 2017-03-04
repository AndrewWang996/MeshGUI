function newVertices = reorient(meshname, vertices)
    [V,F] = getMesh(meshname);
    
    anchorIndices = getAnchorIndices(meshname);
    originalAnchorPositions = V(anchorIndices,1:2);
    anchorPositions = vertices(anchorIndices,1:2);
    
    A = inv(anchorPositions) * originalAnchorPositions;
    newVertices = vertices(:,1:2) * A;

end