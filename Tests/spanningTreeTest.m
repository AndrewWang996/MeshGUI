function spanningTreeTest()
    meshname = 'red_dragon';
    [endNodes, weights, predecessor] = getSpanningTree(meshname);
    display('done');
end