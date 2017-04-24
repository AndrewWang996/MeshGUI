function meshToGraphTest()
    addpath Helpers/OptimizationHelpers;
    [vertices,faces] = getMesh('simple');
    G = meshToGraph(vertices, faces);
    plot(G);
    
end