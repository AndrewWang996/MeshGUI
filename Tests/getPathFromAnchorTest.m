function getPathFromAnchorTest()
    addpath Helpers/OptimizationHelpers;
    meshname = 'red_dragon';
    [V,F] = getMesh(meshname);
    hold on;
    
    
    trimesh(F, V(:,1), V(:,2));
    
    path = getPathFromAnchor(meshname, 2);
    pathVertices = V(path,:);
    
    quiver(pathVertices(1:end-1,1), pathVertices(1:end-1,2), ...
        pathVertices(2:end,1) - pathVertices(1:end-1,1), ...
        pathVertices(2:end,2) - pathVertices(1:end-1,2));

    display(path);
    
    hold off;
    
end