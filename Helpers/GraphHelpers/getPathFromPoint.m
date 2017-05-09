%% returns a list of vertex indices in the mesh
% that starts with the first anchor point in the mesh
% and ends at the given vertex (which is given as an integer index)
function path = getPathFromPoint(meshname, fromIndex, toIndex)
    [vertices, faces] = getMesh(meshname);
    G = meshToGraph(vertices, faces);
    
    path = shortestpath(G, fromIndex, toIndex);
end