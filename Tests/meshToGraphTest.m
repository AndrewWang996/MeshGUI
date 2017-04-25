function meshToGraphTest()
    addpath Helpers/OptimizationHelpers;
    [vertices,faces] = getMesh('vert_bar');
    G = meshToGraph(vertices, faces);
    plot(G);
    
end