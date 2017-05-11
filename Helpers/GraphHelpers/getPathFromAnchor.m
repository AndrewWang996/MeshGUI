%% returns a list of vertex indices in the mesh
% that starts with the first anchor point in the mesh
% and ends at the given vertex (which is given as an integer index)
function path = getPathFromAnchor(meshname, vertexIndex)
    anchorIndices = getAnchorIndices(meshname);
    anchorIndex = anchorIndices(1,:);     % assuming we are using the first anchor
    
    path = getPathFromPoint(anchorIndex, vertexIndex);
end